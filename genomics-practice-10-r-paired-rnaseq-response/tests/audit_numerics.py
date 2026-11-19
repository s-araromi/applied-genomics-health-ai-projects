#!/usr/bin/env python3
"""Optional independent Python audit of the R exports; no data are generated.

Run from the project directory after installing requirements-validation.txt.
The main R analysis and R tests do not require Python or these packages.
"""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import norm, false_discovery_control

ROOT=Path('.')


def rows(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))


def main():
    result=rows('outputs/paired_differential_expression.csv')
    ids=[r['gene_id'] for r in result]
    beta=np.array([float(r['log2_fold_change']) for r in result])
    se=np.array([float(r['lfc_se']) for r in result])
    p0=np.array([float(r['p_point_null']) for r in result])
    p1=np.array([float(r['p_above_1_5_fold']) for r in result])
    q1=np.array([float(r['q_above_1_5_fold']) for r in result])
    expected0=2*norm.sf(abs(beta/se))
    expected1=2*norm.sf(np.maximum((abs(beta)-np.log2(1.5))/se,0))
    # Relative comparison also checks extremely small probabilities; parsing
    # with float avoids rounded long-decimal CSV-reader representations.
    np.testing.assert_allclose(p0,expected0,rtol=1e-10,atol=1e-300)
    np.testing.assert_allclose(p1,expected1,rtol=1e-10,atol=1e-300)
    qref=false_discovery_control(p1,method='bh')
    np.testing.assert_allclose(q1,qref,rtol=1e-10,atol=1e-300)
    source={r['gene_id']:r for r in rows('data/counts.csv')}
    samples=rows('data/samples.csv');runs=[r['run'] for r in samples]
    counts=np.array([[int(source[g][run]) for run in runs] for g in ids])
    log_counts=np.log(counts[(counts>0).all(axis=1)])
    factors=np.exp(np.median(log_counts-log_counts.mean(axis=1)[:,None],axis=0))
    qc=rows('outputs/sample_qc.csv')
    np.testing.assert_allclose(factors,[float(r['size_factor']) for r in qc],rtol=1e-12)
    assert [r['run'] for r in qc]==runs
    selected=q1<.05
    assert selected.sum()==600
    assert (beta[selected]>0).sum()==352
    assert (beta[selected]<0).sum()==248
    print('Independent audit passed: normalization, both Wald tests, full-family BH, and selection counts.')
    print('Maximum absolute discrepancies:')
    for label,actual,reference in [('size_factor',factors,np.array([float(r['size_factor']) for r in qc])),
                                   ('point_p',p0,expected0),('threshold_p',p1,expected1),('threshold_q',q1,qref)]:
        print(f'{label}: {np.max(abs(actual-reference)):.3g}')


if __name__=='__main__':main()
