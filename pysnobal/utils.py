def hrs2min(x): return x * 60


def min2sec(x): return x * 60


C_TO_K = 273.16
FREEZE = C_TO_K


def SEC_TO_HR(x): return x / 3600.0

# Kelvin to Celcius


def K_TO_C(x): return x - FREEZE


def check_range(value, min_val, max_val, descrip):
    """
    Check the range of the value
    Args:
        value: value to check
        min_val: minimum value
        max_val: maximum value
        descrip: short description of input
    Returns:
        True if within range
    """
    if (value < min_val) or (value > max_val):
        raise ValueError("%s (%f) out of range: %f to %f",
                         descrip, value, min_val, max_val)
    return True
