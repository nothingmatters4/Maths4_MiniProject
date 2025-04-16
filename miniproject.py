from sympy import *
from sympy.abc import x, a, z
from sympy.integrals.transforms import inverse_laplace_transform

def residue_theorem_integral():
    print("Enter the function in terms of x (example: 1/(x**2 + 1)):")
    fx_input = input("f(x) = ")
    try:
        f = sympify(fx_input)
    except:
        print("Invalid function format.")
        return
    print("Enter the lower limit of integration (use -oo for -infinity):")
    a_lim = sympify(input("Lower limit = "))
    print("Enter the upper limit of integration (use oo for infinity):")
    b_lim = sympify(input("Upper limit = "))
    try:
        result_sympy = integrate(f, (x, a_lim, b_lim))
        print(f"\n[SymPy] Result of the definite integral: {result_sympy}")
    except:
        print("Could not compute the integral using SymPy.")
    if a_lim == -oo and b_lim == oo:
        z_real = symbols('z')
        fz = f.subs(x, z_real)
        poles = singularities(fz, z_real)
        residues = []
        for p in poles:
            if im(p) > 0:
                res = residue(fz, z_real, p)
                residues.append(res)
        integral_residue = 2 * pi * I * sum(residues)
        print(f"[Residue Theorem] Integral over real axis = {integral_residue.evalf()}")
    elif a_lim == 0 and b_lim == 2*pi:
        if f.has(cos(x)):
            try:
                num, denom = fraction(f)
                if num.has(cos(2*x)) and denom.has(cos(x)):
                    print("\n[Residue Method - Unit Circle] Attempting complex substitution...")
                    cos_theta = (z + 1/z)/2
                    cos_2theta = (z**2 + 1/z**2)/2
                    dtheta = 1 / (I * z)
                    integrand_z = (cos_2theta / (1 - 2*a*cos_theta + a**2)) * dtheta
                    integrand_z = simplify(integrand_z)
                    print("Transformed integrand over unit circle |z| = 1:")
                    pprint(integrand_z)
                    poles = singularities(integrand_z, z)
                    inside_poles = [p for p in poles if abs(p.evalf(subs={a: 0.5})) < 1.01]
                    print("\nPoles inside unit circle:")
                    pprint(inside_poles)
                    residues = [residue(integrand_z, z, p) for p in inside_poles]
                    for p, r in zip(inside_poles, residues):
                        print(f"Residue at z = {p} → {r}")
                    result = 2 * pi *I* sum(residues)
                    print(f"\n[Final Result] ∫₀²π cos(2θ)/(1 - 2a cos θ + a²) dθ = {result.simplify()}")
                else:
                    print("\nFunction is not of the form cos(2x)/(1 - 2a cos x + a^2); skipping unit circle method.")
            except Exception as e:
                print(f"\nCould not apply unit circle residue method: {e}")
        else:
            print("\nFunction does not involve cos(x); unit circle residue method skipped.")
    else:
        print("\nResidue method demonstration skipped (only applies to ∫ from -∞ to ∞ or [0, 2π] with trigonometric forms).")

residue_theorem_integral()
