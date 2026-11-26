"""Composable regression baselines for simulation experiments.

NumPy is required; scikit-learn and PyTorch are optional training dependencies.
This package never replaces a simulator force or transition law implicitly.
"""
from .data import Dataset, DataSplit, RegressionTask, Standardizer, grouped_split
from .methods import available_methods, register_method
from .efficiency import EfficiencyReport, evaluate_learning_efficiency
from .platform import (ARTIFACT_SCHEMA_VERSION, TrainedArtifact, compare,
                       load_artifact, regression_metrics, train)

__all__ = ["Dataset", "DataSplit", "RegressionTask", "Standardizer", "grouped_split",
           "available_methods", "register_method", "TrainedArtifact", "train", "compare",
           "load_artifact", "regression_metrics", "ARTIFACT_SCHEMA_VERSION",
           "EfficiencyReport", "evaluate_learning_efficiency"]
