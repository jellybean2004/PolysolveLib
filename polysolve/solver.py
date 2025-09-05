from cmath import sqrt
import math

from cowsay import cow

CBRT_UNITY_IM = sqrt(3)/2 * 1j

def quadratic(a: float, b: float, c: float) -> tuple[complex, complex]:
    """
    Solves the roots of a quadratic equation
    
    Uses the quadratic formula. Results must be real
    
    Parameters
    ----------
    a
        :math: `x^2` coefficient
    b
        :math: `x` coefficient
    c
        Constant coefficient
    
    Returns
    -------
    tuple[complex, complex]
        Quadratic roots
        
    Examples
    --------
    >>> quadratic(3., 0., -1.)
    ((0.5773502691896257+0j), (-0.5773502691896257+0j))
    >>> quadratic(1, 2, 0)
    (0j, (-2+0j))
    >>> quadratic(1, 11, 28)
    ((-4+0j), (-7+0j))
    
    Raises
    ------
    ValueError
        Discriminant < 0 implying imaginary root
        
    Notes
    -----
    The quadratic formula is given by:
    
    .. math::
        
        x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
        
    See Also
    --------
    cubic : Solves the roots of a cubic equation
    
    References
    ----------
    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
           expert systems ...
    """
    det = b**2 - (4*a*c)

    if math.isclose(det, 0):
        cow("Degenerate MOOoo-ts")

    return ((-b + sqrt(det)) / (2*a), (-b - sqrt(det)) / (2*a))

def cubic(a: float, b: float, c: float, d: float) -> tuple[complex, complex, complex]:
    """
    Solves the roots of a cubic equation
    
    Uses the cubic formula. Results must be real
    
    Parameters
    ----------
    a : float
        :math: `x**3` coefficient
    b : float
        :math: `x**2` coefficient
    c : float
        :math: `x` coefficient
    d : float
        Constant coefficient
        
    Returns
    -------
    tuple[complex, complex, complex]
        Cubic roots
    
    Raises
    ------
    ValueError
        If the cubic equation has no real roots
        
    Notes
    -----
    The cubic formula is given by:
    
    .. math::
        
        x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
        
    See Also
    --------
    quadratic : Solves the roots of a quadratic equation
    
    References
    ----------
    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
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