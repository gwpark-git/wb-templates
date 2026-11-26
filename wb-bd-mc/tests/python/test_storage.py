import importlib.util
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import threading

import pytest
from bdmc.storage import DataBudget,DataBudgetExceeded

spec=importlib.util.spec_from_file_location('bdmc_storage',Path(__file__).resolve().parents[2]/'tools/check_storage.py')
storage=importlib.util.module_from_spec(spec)
spec.loader.exec_module(storage)


def test_storage_stops_before_projected_limit(tmp_path):
    (tmp_path/'source.txt').write_text('retained source')
    report=storage.check([tmp_path],planned_bytes=5_000_000_000)
    assert report['status']=='stop_before_expansion'
    assert (tmp_path/'source.txt').read_text()=='retained source'
    assert storage.check([tmp_path])['status']=='ok'


def test_storage_nested_roots_rejected(tmp_path):
    child=tmp_path/'child';child.mkdir()
    with pytest.raises(ValueError,match='nonoverlapping'):
        storage.check([tmp_path,child])


def test_additional_budget_excludes_existing_files_and_refuses_reset(tmp_path):
    data=tmp_path/'runs';data.mkdir()
    (data/'existing.bin').write_bytes(b'x'*100_000)
    baseline=tmp_path/'baseline.json'
    storage.create_data_baseline([data],baseline)
    assert storage.check_additional(baseline)['additional_bytes']==0
    assert storage.check_additional(baseline,planned_bytes=8_000_000_000)['status']=='ok'
    assert storage.check_additional(baseline,planned_bytes=9_000_000_000)['status']=='stop_before_expansion'
    (data/'new.bin').write_bytes(b'y'*25_000)
    assert storage.check_additional(baseline)['additional_bytes']>=25_000
    with pytest.raises(FileExistsError):
        storage.create_data_baseline([data],baseline)


def test_removing_old_data_cannot_mask_another_root_growth(tmp_path):
    a,b=tmp_path/'a',tmp_path/'b';a.mkdir();b.mkdir()
    (a/'old').write_bytes(b'x'*100_000)
    baseline=tmp_path/'baseline.json';storage.create_data_baseline([a,b],baseline)
    (a/'old').unlink();(b/'new').write_bytes(b'y'*25_000)
    assert storage.check_additional(baseline)['additional_bytes']>=25_000


def test_budget_wrapper_validates_output_scope_and_planned_generation(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'baseline.json';storage.create_data_baseline([root],baseline)
    budget=DataBudget(baseline)
    assert budget.require_output(root/'new',planned_bytes=1000)['status']=='ok'
    with pytest.raises(ValueError,match='outside'):
        budget.require_output(tmp_path/'outside',planned_bytes=1000)
    with pytest.raises(DataBudgetExceeded):
        budget.require_output(root/'new',planned_bytes=10_000_000_000)
    assert not (root/'new').exists()


def test_concurrent_writers_recheck_under_one_shared_lock(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';storage.create_data_baseline([root],baseline)
    record=json.loads(baseline.read_text())
    record.update(soft_limit_bytes=16_384,additional_limit_bytes=20_480)
    baseline.write_text(json.dumps(record))
    start=threading.Barrier(2)
    def write(index):
        budget=DataBudget(baseline)
        start.wait(timeout=5)
        try:
            with budget.writer(root/f'data{index}',planned_bytes=12_288):
                (root/f'data{index}').write_bytes(b'x'*12_288)
            return 'written'
        except DataBudgetExceeded:
            return 'stopped'
    with ThreadPoolExecutor(max_workers=2) as pool:
        result=list(pool.map(write,[0,1]))
    assert sorted(result)==['stopped','written']
    assert DataBudget(baseline).check()['additional_bytes']==12_288


def test_nested_writers_share_lock_and_exception_releases_it(tmp_path):
    root=tmp_path/'runs';root.mkdir()
    baseline=tmp_path/'budget.json';storage.create_data_baseline([root],baseline)
    budget=DataBudget(baseline)
    with pytest.raises(RuntimeError):
        with budget.writer(root/'a',planned_bytes=1000):
            with DataBudget(baseline).writer(root/'b',planned_bytes=1000):
                raise RuntimeError('simulated write error')
    with budget.writer(root/'c',planned_bytes=1000):
        (root/'c').write_bytes(b'ok')
