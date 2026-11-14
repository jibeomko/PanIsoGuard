# SQANTI-SIM's pre-trained NanoSim models were pickled with scikit-learn 1.2, whose distance class
# scikit-learn >= 1.3 renamed (EuclideanDistance -> EuclideanDistance64). Put this directory on
# PYTHONPATH for the NanoSim step and the pickles load; the KernelDensity objects sample as before.
try:
    import sklearn.metrics._dist_metrics as _m
    if not hasattr(_m, "EuclideanDistance") and hasattr(_m, "EuclideanDistance64"):
        _m.EuclideanDistance = _m.EuclideanDistance64
except Exception:
    pass
