"""Small regression methods with numeric-only inference state.

New methods register trusted fit/predict callables. They must keep all inference
parameters in finite numeric arrays and a JSON configuration; no pickle loader
is part of this contract.
"""
from __future__ import annotations

from dataclasses import dataclass
import importlib.util
from typing import Any, Callable

import numpy as np


@dataclass(frozen=True)
class Method:
    fit: Callable[..., tuple[dict[str, np.ndarray], dict[str, Any]]]
    predict: Callable[..., tuple[np.ndarray, np.ndarray | None]]
    dependency: str | None = None
    uncertainty: str | None = None
    predict_std: Callable | None = None


_METHODS: dict[str, Method] = {}


def register_method(name: str, *, fit: Callable, predict: Callable,
                    dependency: str | None = None, uncertainty: str | None = None,
                    predict_std: Callable | None = None) -> None:
    """Register an explicit local extension; existing names cannot be replaced."""
    # Registration executes trusted caller-provided code, not code deserialized
    # from a model file. Keep method choice explicit and retain the numeric-state
    # contract so saved inference need not restore a training-library object.
    if not isinstance(name, str) or not name or name in _METHODS:
        raise ValueError("method name must be nonempty and not already registered")
    if not callable(fit) or not callable(predict):
        raise ValueError("fit and predict must be callable")
    if predict_std is not None and (not callable(predict_std) or not uncertainty):
        raise ValueError("predict_std must be callable and declare an uncertainty interpretation")
    _METHODS[name] = Method(fit, predict, dependency, uncertainty, predict_std)


def available_methods() -> dict[str, dict[str, Any]]:
    return {name: {"available": method.dependency is None
                  or importlib.util.find_spec(method.dependency) is not None,
                   "dependency": method.dependency, "uncertainty": method.uncertainty}
            for name, method in _METHODS.items()}


def get_method(name: str, *, training: bool = False) -> Method:
    if name not in _METHODS:
        raise ValueError(f"unknown method {name!r}; registered methods: {list(_METHODS)}")
    method = _METHODS[name]
    if training and method.dependency and importlib.util.find_spec(method.dependency) is None:
        raise ImportError(f"training {name} requires optional dependency {method.dependency}")
    return method


def _check_keys(config: dict[str, Any], allowed: set[str]) -> None:
    unknown = set(config) - allowed
    if unknown:
        raise ValueError(f"unknown configuration keys: {sorted(unknown)}")


def _positive(value: Any, name: str, *, allow_zero: bool = False) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0 or (result == 0 and not allow_zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return result


def _integer(value: Any, name: str, *, maximum: int | None = None) -> int:
    result = int(value)
    if result != value or isinstance(value, bool) or result <= 0 or (maximum and result > maximum):
        raise ValueError(f"{name} must be a positive integer" + (f" <= {maximum}" if maximum else ""))
    return result


def _ridge_fit(X: np.ndarray, y: np.ndarray, config: dict, seed: int):
    _check_keys(config, {"alpha"})
    alpha = _positive(config.get("alpha", 1e-6), "alpha", allow_zero=True)
    design = np.column_stack((np.ones(len(X)), X))
    # Augmented least squares also handles rank-deficient designs and alpha=0.
    penalty = np.diag(np.r_[0.0, np.full(X.shape[1], np.sqrt(alpha))])
    coefficients = np.linalg.lstsq(np.vstack((design, penalty)),
                                   np.vstack((y, np.zeros((X.shape[1] + 1, y.shape[1])))),
                                   rcond=None)[0]
    return {"coefficients": coefficients}, {"alpha": alpha}


def _ridge_predict(X: np.ndarray, state: dict, config: dict):
    coefficients = state["coefficients"]
    if coefficients.ndim != 2 or coefficients.shape[0] != X.shape[1] + 1:
        raise ValueError("invalid ridge coefficient shape")
    return np.column_stack((np.ones(len(X)), X)) @ coefficients, None


def _gp_fit(X: np.ndarray, y: np.ndarray, config: dict, seed: int):
    _check_keys(config, {"length_scale", "signal_variance", "noise_variance"})
    length_scale = _positive(config.get("length_scale", 1.0), "length_scale")
    signal_variance = _positive(config.get("signal_variance", 1.0), "signal_variance")
    noise_variance = _positive(config.get("noise_variance", 1e-4), "noise_variance")
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, RBF

    kernel = ConstantKernel(signal_variance, constant_value_bounds="fixed") * RBF(
        length_scale, length_scale_bounds="fixed")
    model = GaussianProcessRegressor(kernel=kernel, alpha=noise_variance,
                                     optimizer=None, normalize_y=False, random_state=seed)
    model.fit(X, y)
    alpha = model.alpha_
    if alpha.ndim == 1:
        alpha = alpha[:, None]
    return {"X_train": X.copy(), "alpha": alpha, "cholesky": model.L_}, {
        "length_scale": length_scale, "signal_variance": signal_variance,
        "noise_variance": noise_variance}


def _gp_cross_kernel(X: np.ndarray, state: dict, config: dict):
    X_train = state["X_train"]
    alpha, cholesky = state["alpha"], state["cholesky"]
    n = len(X_train)
    if (X_train.ndim != 2 or X_train.shape[1] != X.shape[1] or alpha.ndim != 2
            or len(alpha) != n or cholesky.shape != (n, n)):
        raise ValueError("invalid Gaussian process state shapes")
    squared_distance = ((X[:, None, :] - X_train[None, :, :]) ** 2).sum(axis=2)
    kernel = config["signal_variance"] * np.exp(-0.5 * squared_distance / config["length_scale"] ** 2)
    return kernel


def _gp_predict(X: np.ndarray, state: dict, config: dict):
    kernel = _gp_cross_kernel(X, state, config)
    return kernel @ state["alpha"], None


def _gp_std(X: np.ndarray, state: dict, config: dict):
    # This is conditional latent-function spread under the chosen kernel. It
    # does not estimate between-trajectory physical noise, timestep bias or
    # calibrated coverage of a new scientific regime.
    kernel = _gp_cross_kernel(X, state, config)
    cholesky = state["cholesky"]
    if np.any(np.diag(cholesky) <= 0):
        raise ValueError("Gaussian process Cholesky diagonal must be positive")
    # Forward substitution uses the stored factor in O(n_train**2 * n_query).
    # A generic np.linalg.solve would refactor this triangular matrix cubically.
    projection = np.empty_like(kernel.T)
    for i in range(len(cholesky)):
        projection[i] = (kernel[:, i] - cholesky[i, :i] @ projection[:i]) / cholesky[i, i]
    variance = np.maximum(0.0, config["signal_variance"] - (projection ** 2).sum(axis=0))
    return np.repeat(np.sqrt(variance)[:, None], state["alpha"].shape[1], axis=1)


def _mlp_fit(X: np.ndarray, y: np.ndarray, config: dict, seed: int):
    _check_keys(config, {"hidden_layers", "epochs", "learning_rate", "weight_decay"})
    hidden = tuple(_integer(width, "hidden width", maximum=4096)
                   for width in config.get("hidden_layers", (32, 32)))
    if not hidden or len(hidden) > 8:
        raise ValueError("hidden_layers must have one to eight widths")
    epochs = _integer(config.get("epochs", 250), "epochs")
    learning_rate = _positive(config.get("learning_rate", 0.01), "learning_rate")
    weight_decay = _positive(config.get("weight_decay", 0.0), "weight_decay", allow_zero=True)
    import torch

    # CPU-only prototype. fork_rng restores the caller's CPU RNG state.
    with torch.random.fork_rng(devices=[]):
        torch.manual_seed(seed)
        widths = (X.shape[1], *hidden, y.shape[1])
        layers = []
        for i in range(len(widths) - 1):
            layers.append(torch.nn.Linear(widths[i], widths[i + 1], dtype=torch.float64))
            if i < len(widths) - 2:
                layers.append(torch.nn.Tanh())
        model = torch.nn.Sequential(*layers)
        optimiser = torch.optim.Adam(model.parameters(), lr=learning_rate,
                                     weight_decay=weight_decay)
        input_tensor = torch.tensor(X, dtype=torch.float64)
        target_tensor = torch.tensor(y, dtype=torch.float64)
        for _ in range(epochs):
            optimiser.zero_grad()
            loss = torch.nn.functional.mse_loss(model(input_tensor), target_tensor)
            if not bool(torch.isfinite(loss)):
                raise ValueError("MLP loss became nonfinite")
            loss.backward()
            optimiser.step()
        state = {}
        i = 0
        for layer in model:
            if isinstance(layer, torch.nn.Linear):
                state[f"weight_{i}"] = layer.weight.detach().numpy().copy()
                state[f"bias_{i}"] = layer.bias.detach().numpy().copy()
                i += 1
    return state, {"hidden_layers": list(hidden), "epochs": epochs,
                   "learning_rate": learning_rate, "weight_decay": weight_decay,
                   "activation": "tanh", "device": "cpu", "dtype": "float64"}


def _mlp_predict(X: np.ndarray, state: dict, config: dict):
    # Exported weights support plain NumPy inference after PyTorch training.
    # Keeping that boundary explicit tests the saved model actually used for
    # prediction, rather than relying on a live in-memory training object.
    result = X
    n_layers = len(config["hidden_layers"]) + 1
    for i in range(n_layers):
        weight, bias = state[f"weight_{i}"], state[f"bias_{i}"]
        if (weight.ndim != 2 or bias.shape != (weight.shape[0],)
                or weight.shape[1] != result.shape[1]):
            raise ValueError("invalid MLP parameter shapes")
        result = result @ weight.T + bias
        if i < n_layers - 1:
            result = np.tanh(result)
    return result, None


register_method("ridge", fit=_ridge_fit, predict=_ridge_predict)
register_method("gaussian_process", fit=_gp_fit, predict=_gp_predict, dependency="sklearn",
                uncertainty="latent_function_standard_deviation_uncalibrated", predict_std=_gp_std)
register_method("mlp", fit=_mlp_fit, predict=_mlp_predict, dependency="torch")
