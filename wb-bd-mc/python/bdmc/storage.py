"""Scoped additional-data accounting for unattended BD–MC generation.

The fixed baseline excludes existing files and environments. Checks are
preflights, not a filesystem quota; callers must bound and recheck writes.
"""
from datetime import datetime,timezone
from contextlib import contextmanager
import json
from pathlib import Path
import subprocess
import threading

_WRITER_LOCKS = {}
_WRITER_LOCKS_GUARD = threading.Lock()
_WRITER_LOCAL = threading.local()

def allocated_bytes(path):
    result=subprocess.run(['du','-skx',str(path)],capture_output=True,text=True,check=True)
    return int(result.stdout.split()[0])*1024


def check(roots,planned_bytes=0,soft_limit=4_500_000_000,hard_limit=5_000_000_000):
    """Legacy total-usage diagnostic; no longer the researcher's active policy."""
    roots=[Path(root).resolve(strict=True) for root in roots]
    if planned_bytes<0 or not 0<soft_limit<hard_limit:
        raise ValueError('nonnegative planned bytes and 0<soft<hard limits required')
    # Nested roots would double-count; reject rather than silently changing scope.
    for i,root in enumerate(roots):
        if any(root==other or root in other.parents or other in root.parents for other in roots[:i]):
            raise ValueError('storage roots must be distinct and nonoverlapping')
    usage={str(root):allocated_bytes(root) for root in roots}
    projected=sum(usage.values())+planned_bytes
    return {'created_at':datetime.now(timezone.utc).isoformat(),'allocated_bytes':usage,
            'total_bytes':sum(usage.values()),'planned_bytes':planned_bytes,'projected_bytes':projected,
            'soft_limit_bytes':soft_limit,'hard_limit_bytes':hard_limit,
            'status':'ok' if projected<soft_limit else 'stop_before_expansion',
            'measurement':'du allocated bytes, scoped paths only; no unrelated workspace traversal'}


def create_data_baseline(roots,path):
    """Never silently reset the starting point and grant another data allowance."""
    path=Path(path)
    if path.exists():
        raise FileExistsError(path)
    usage=check(roots)['allocated_bytes']
    record={'schema':'bdmc.additional-data-budget.v1',
            'created_at':datetime.now(timezone.utc).isoformat(),
            'baseline_allocated_bytes':usage,'additional_limit_bytes':10_000_000_000,
            'soft_limit_bytes':9_000_000_000,
            'scope':'additional generated on-disk data from this baseline; existing data/environments excluded',
            'measurement':'scoped data roots, allocated KiB rounded to bytes; no filesystem quota'}
    with path.open('x') as stream:
        stream.write(json.dumps(record,indent=2)+'\n')
    return record


def check_additional(baseline_path,planned_bytes=0):
    if planned_bytes<0:
        raise ValueError('planned bytes must be nonnegative')
    baseline=json.loads(Path(baseline_path).read_text())
    if baseline['schema']!='bdmc.additional-data-budget.v1':
        raise ValueError('unsupported data budget schema')
    before=baseline['baseline_allocated_bytes']
    current={root:allocated_bytes(Path(root)) if Path(root).exists() else 0 for root in before}
    # Removing baseline data from one root must not mask growth in another.
    growth={root:max(0,current[root]-size) for root,size in before.items()}
    projected=sum(growth.values())+planned_bytes
    return {'created_at':datetime.now(timezone.utc).isoformat(),
            'budget_type':'additional_generated_data','baseline':str(Path(baseline_path).resolve()),
            'baseline_created_at':baseline['created_at'],'baseline_allocated_bytes':before,
            'current_allocated_bytes':current,'additional_bytes_by_root':growth,
            'additional_bytes':sum(growth.values()),'planned_bytes':planned_bytes,
            'projected_additional_bytes':projected,'soft_limit_bytes':baseline['soft_limit_bytes'],
            'hard_limit_bytes':baseline['additional_limit_bytes'],
            'status':'ok' if projected<baseline['soft_limit_bytes'] else 'stop_before_expansion',
            'measurement':baseline['measurement']}


class DataBudgetExceeded(RuntimeError):
    pass


class DataBudget:
    def __init__(self, baseline):
        self.path=Path(baseline).resolve(strict=True)
        # Validate by reading the current status without resetting anything.
        self.check()

    def check(self, planned_bytes=0):
        return check_additional(self.path,planned_bytes)

    def require_output(self, output_path, planned_bytes=0):
        output=Path(output_path).resolve()
        report=self.check(planned_bytes)
        roots=[Path(root) for root in report['baseline_allocated_bytes']]
        if not any(output==root or root in output.parents for root in roots):
            raise ValueError('output is outside the fixed data-budget roots; register its scope before generating')
        if report['status']!='ok':
            raise DataBudgetExceeded(
                f"additional data plus planned output reaches the {report['soft_limit_bytes']}-byte stop point")
        return report

    @contextmanager
    def writer(self, output_path, planned_bytes=0):
        """Serialize cooperating disk writers and recheck within the shared lock.

        Compute remains parallel. The lock is released by the OS on process
        exit; no stale reservation can grant another allowance. Same-thread
        nested dataset/checkpoint writes use one lock for the baseline.
        """
        # A preflight outside the lock races with another authorized writer.
        # Thread and file locks share one baseline key; rechecking after entry
        # prevents cooperating processes from spending the same remaining space.
        # Writers bypassing this API are outside the guarantee, not a disk quota.
        import fcntl
        key=str(self.path)
        with _WRITER_LOCKS_GUARD:
            lock=_WRITER_LOCKS.setdefault(key,threading.RLock())
        with lock:
            held=getattr(_WRITER_LOCAL,'held',None)
            if held is None:
                held={};_WRITER_LOCAL.held=held
            if key in held:
                yield self.require_output(output_path,planned_bytes)
                return
            lock_path=self.path.with_name(self.path.name+'.write.lock')
            with lock_path.open('a') as handle:
                fcntl.flock(handle.fileno(),fcntl.LOCK_EX)
                held[key]=handle
                try:
                    yield self.require_output(output_path,planned_bytes)
                finally:
                    held.pop(key,None)
                    fcntl.flock(handle.fileno(),fcntl.LOCK_UN)
