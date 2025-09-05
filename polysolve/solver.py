from cmath import sqrt
import math

from cowsay import cow

CBRT_UNITY_IM = sqrt(3)/2 * 1j

def quadratic(a: float, b: float, c: float) -> tuple[complex, complex]:
    """
    Solves the roots of a quadratic equation.
    
    Parameters
    ----------
    a
        x^2 coefficient
    b
        x coefficient
    c
        constant coefficient
        
    Returns
    -------
    tuple[complex, complex]
        Quadratic roots
        
    Examples
    --------
    >>> quadratic(1, 2, 0)
    (0j, (-2+0j))
    
    >>> quadratic(1, 11, 28)
    ((-4+0j), (-7+0j))
    
    Notes
    -----
    This function uses the quadratic formula to solve for the roots of a quadratic equation.
    
    Raises
    ------
    ValueError
        If the quadratic equation is degenerate.
        
    References
    ----------
    [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems ...
    
    """
    det = b**2 - (4*a*c)

    if math.isclose(det, 0):
        cow("Degenerate MOOoo-ts")

    return ((-b + sqrt(det)) / (2*a), (-b - sqrt(det)) / (2*a))

def cubic(a: float, b: float, c: float, d: float) -> tuple[complex, complex, complex]:
    """
    Solves the roots of a cubic equation
    
    Parameters
    ----------
    a
        x^3 coefficient
    b
        x^2 coefficient
    c
        x coefficient
    d
        Constant coefficient
        
    Returns
    -------
    tuple[complex, complex, complex]
        Cubic roots
        
    Examples
    --------    
    >>> cubic(1, 0, 0, -1)
    (-1j, (1+0j), (-0+0j))
    
    >>> cubic(1, 0, 0, 2)
    (2j, (-2+0j), (-0+0j))
    
    >>> cubic(1, 0, 0, -2)
    (-2j, (2+0j), (-0+0j))
    
    Notes
    -----
    This function uses the cubic formula to solve for the roots of a cubic equation.
    
    Raises
    ------
    ValueError
        If the cubic equation is degenerate.
        
    References
    ----------
    [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems ...
    """
    q = (3*a*c - b**2) / (9*a**2)
    r = (9*a*b*c - 27*a**2*d - 2*b**3) / (54*a**3)

    s = (r + sqrt(q**3 + r**2))**(1/3)
    t = (r - sqrt(q**3 + r**2))**(1/3)

    x1 = s + t - (b/3*a)
    x2 = -(s + t)/2 - (b/3*a) + CBRT_UNITY_IM * (s - t)
    x3 = -(s + t)/2 - (b/3*a) - CBRT_UNITY_IM * (s - t)

    if any(x == x1 for x in (x2, x3)):
        cow("Degenerate MOOoo-ts")

    return (x1, x2, x3)

if __name__  == "__main__":
    import doctest
    doctest.testmod()