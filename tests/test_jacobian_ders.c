#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/logging.h>

#include "types.h"


/*
   Forward declarations for the internal stencil helpers in
   source/jacobian_ders.c. These functions implement individual
   stencil blocks of the elliptic Newton Jacobian and are called
   by the dispatcher ell_interior_point_ders. They are file-scope
   externals (not static), but are NOT exposed in
   include/jacobian_ders.h. Declared here so the unit tests can
   call them directly without going through the dispatcher.

   The three tested below are the three shape classes:
     - im1_jm1: diagonal corner block (representative of the four
                diagonal corners im1_jm1, im1_jp1, ip1_jm1, ip1_jp1
                which differ only in sign).
     - i_j:    center block (diagonal of the system).
     - ip1_j:  cardinal block (representative of the four cardinals
                im1_j, ip1_j, i_jm1, i_jp1 which differ in sign and
                axis only).
*/
ell_jacobian_2D ell_interior_im1_jm1 (long double, long double, coeffs_1);
ell_jacobian_2D ell_interior_i_j (long double, long double, coeffs_1);
ell_jacobian_2D ell_interior_ip1_j (long double, grid_der_2D, coeffs_1, grid_dder_2D);




/*
    Test ell_interior_im1_jm1: derivative of F(u(i, j)) with respect
    to the diagonal corner neighbor (i - 1, j - 1).

    The diagonal corner blocks come from the central-FD discretization
    of the cross derivative r_xieta in the elliptic operator
    F = alpha r_xixi - 2 beta r_xieta + gamma r_etaeta. The four-point
    cross stencil places (i - 1, j - 1) with weight +1/(4 dxi deta),
    so its contribution to F(i, j) is -2 beta * 1/(4 dxi deta) =
    -beta/(2 dxi deta), and the derivative wrt the neighbor reproduces
    that coefficient on the diagonal of the 2x2 block. Cross
    components (ddx.y, ddy.x) are zero because the F_x equation has
    no dependence on y_{i-1,j-1} (and vice versa) in this stencil.

    Closed form:
        ders.ddx.x = -beta/(2 dxi deta)        (= dF_x/dx_{i-1,j-1})
        ders.ddy.y = -beta/(2 dxi deta)        (= dF_y/dy_{i-1,j-1})
        ders.ddx.y = 0
        ders.ddy.x = 0

    Concrete inputs:
        dxi  = 0.5
        deta = 0.25
        coeff_ij.beta  = 0.5
        (alpha, gamma, J set to filler values to confirm they are
         not read by this stencil)
    Expected diagonal: -0.5/(2 * 0.5 * 0.25) = -2.0
*/
Test(jacobian_ders, ell_interior_im1_jm1)
{
    // Inputs
    long double             dxi = 0.5L;
    long double             deta = 0.25L;
    coeffs_1                coeff_ij;

    // Output and expected
    ell_jacobian_2D         ders;
    long double             expected_diag;


    // Populate coeffs: only beta is read; the other fields are
    // set to filler values that should not appear in the output
    coeff_ij.alpha = 2.0L;
    coeff_ij.beta = 0.5L;
    coeff_ij.gamma = 3.0L;
    coeff_ij.J = 1.0L;

    // Closed-form expected value
    expected_diag = -coeff_ij.beta/(2.0L * dxi * deta);


    // Call function
    ders = ell_interior_im1_jm1 (dxi, deta, coeff_ij);


    // Diagonal entries match closed form
    cr_assert_float_eq (ders.ddx.x, expected_diag, 1e-15L,
                        "ddx.x mismatch: got %.15Le, expected %.15Le",
                        ders.ddx.x, expected_diag);
    cr_assert_float_eq (ders.ddy.y, expected_diag, 1e-15L,
                        "ddy.y mismatch: got %.15Le, expected %.15Le",
                        ders.ddy.y, expected_diag);

    // Off-diagonal entries are zero (no F_x, F_y cross-coupling)
    cr_assert_float_eq (ders.ddx.y, 0.0L, 1e-15L,
                        "ddx.y should be 0 for the diagonal-corner block");
    cr_assert_float_eq (ders.ddy.x, 0.0L, 1e-15L,
                        "ddy.x should be 0 for the diagonal-corner block");
}




/*
    Test ell_interior_i_j: derivative of F(u(i, j)) with respect
    to itself (the center of the 9-point stencil).

    The center block comes from the diagonal contributions of the
    two pure second-derivative stencils r_xixi and r_etaeta in the
    elliptic operator. The 3-point central FD for r_xixi places
    (i, j) with weight -2/dxi^2, giving an alpha r_xixi
    contribution of -2 alpha/dxi^2 to F at (i, j). Similarly
    -2 gamma/deta^2 from the r_etaeta term. The cross derivative
    r_xieta has zero weight at the center. So the diagonal entries
    of the 2x2 block sum the two -2 terms and the off-diagonal
    entries are zero (no x/y coupling at the center).

    Closed form:
        ders.ddx.x = -2 (alpha/dxi^2 + gamma/deta^2)
        ders.ddy.y = -2 (alpha/dxi^2 + gamma/deta^2)
        ders.ddx.y = 0
        ders.ddy.x = 0

    Concrete inputs:
        dxi = 0.5,  deta = 0.25
        coeff_ij.alpha = 2.0
        coeff_ij.gamma = 3.0
        (beta, J set to filler values)
    Expected diagonal: -2 * (2/0.25 + 3/0.0625) = -2 * (8 + 48) = -112.0
*/
Test(jacobian_ders, ell_interior_i_j)
{
    // Inputs
    long double             dxi = 0.5L;
    long double             deta = 0.25L;
    coeffs_1                coeff_ij;

    // Output and expected
    ell_jacobian_2D         ders;
    long double             expected_diag;


    // Populate coeffs: only alpha, gamma are read; beta and J are
    // set to filler values that should not appear in the output
    coeff_ij.alpha = 2.0L;
    coeff_ij.beta = 0.5L;
    coeff_ij.gamma = 3.0L;
    coeff_ij.J = 1.0L;

    // Closed-form expected value
    expected_diag = -2.0L * ((coeff_ij.alpha/(dxi * dxi)) +
                             (coeff_ij.gamma/(deta * deta)));


    // Call function
    ders = ell_interior_i_j (dxi, deta, coeff_ij);


    // Diagonal entries match closed form
    cr_assert_float_eq (ders.ddx.x, expected_diag, 1e-15L,
                        "ddx.x mismatch: got %.15Le, expected %.15Le",
                        ders.ddx.x, expected_diag);
    cr_assert_float_eq (ders.ddy.y, expected_diag, 1e-15L,
                        "ddy.y mismatch: got %.15Le, expected %.15Le",
                        ders.ddy.y, expected_diag);

    // Off-diagonal entries are zero (F_x and F_y decouple at the center)
    cr_assert_float_eq (ders.ddx.y, 0.0L, 1e-15L,
                        "ddx.y should be 0 (F_x and F_y decouple at center)");
    cr_assert_float_eq (ders.ddy.x, 0.0L, 1e-15L,
                        "ddy.x should be 0 (F_x and F_y decouple at center)");
}




/*
    Test ell_interior_ip1_j: derivative of F(u(i, j)) with respect
    to the cardinal neighbor (i + 1, j). This is the only block
    shape in the 9-stencil that exercises the full coefficient-
    derivative machinery: the (i + 1, j) neighbor enters F through

      (a) the direct r_xixi stencil weight +alpha/dxi^2, and
      (b) the dependence of alpha, beta, gamma on first
          derivatives, which are themselves FD'd using (i + 1, j).

    Differentiating F = alpha r_xixi - 2 beta r_xieta + gamma r_etaeta
    with respect to x_{i+1,j} (only beta, gamma depend on x_{i+1,j};
    alpha and the second-derivative stencils r_xixi/r_xieta are
    excluded for this neighbor's effect on those terms because
    they don't involve x_{i+1,j}) yields:

      d/dx_{i+1,j} [F_x] = alpha/dxi^2
                          + (d gamma/dx_{i+1,j}) * x_etaeta
                          + (d beta /dx_{i+1,j}) * (-2 x_xieta)

    With d gamma/dx_{i+1,j} = x_xi/dxi  and  d beta/dx_{i+1,j} = x_eta/(2 dxi):

      d/dx_{i+1,j} F_x = (1/dxi) * [alpha/dxi + x_xi*x_etaeta - x_eta*x_xieta]

    Analogous expressions for the other three components follow by
    swapping x <-> y in the appropriate slots.

    Closed form (with mult = 1/dxi):
        ders.ddx.x = mult [ alpha/dxi + x_xi  * x_etaeta - x_eta * x_xieta ]
        ders.ddy.x = mult [             y_xi  * x_etaeta - y_eta * x_xieta ]
        ders.ddx.y = mult [             x_xi  * y_etaeta - x_eta * y_xieta ]
        ders.ddy.y = mult [ alpha/dxi + y_xi  * y_etaeta - y_eta * y_xieta ]

    Concrete inputs (chosen for clean integer arithmetic):
        dxi = 0.5
        coeff_ij.alpha = 2.0   (beta, gamma, J filler)
        x_xi = 1, y_xi = 2, x_eta = 3, y_eta = 4
        x_xixi = 5, y_xixi = 6   (unused for this stencil)
        x_xieta = 7, y_xieta = 8
        x_etaeta = 9, y_etaeta = 10
    With mult = 1/0.5 = 2 and alpha/dxi = 4:
        ders.ddx.x = 2 (4 + 1*9  - 3*7)  = 2 (-8)  = -16
        ders.ddy.x = 2 (    2*9  - 4*7)  = 2 (-10) = -20
        ders.ddx.y = 2 (    1*10 - 3*8)  = 2 (-14) = -28
        ders.ddy.y = 2 (4 + 2*10 - 4*8)  = 2 (-8)  = -16
*/
Test(jacobian_ders, ell_interior_ip1_j)
{
    // Inputs
    long double             dxi = 0.5L;
    coeffs_1                coeff_ij;
    grid_der_2D             dgrid_ij;
    grid_dder_2D            d2grid_ij;

    // Output and expected
    ell_jacobian_2D         ders;
    long double             expected_ddx_x, expected_ddy_x;
    long double             expected_ddx_y, expected_ddy_y;
    long double             mult;


    // coeffs: only alpha is read; beta, gamma, J are fillers
    coeff_ij.alpha = 2.0L;
    coeff_ij.beta = 0.5L;
    coeff_ij.gamma = 3.0L;
    coeff_ij.J = 1.0L;

    // First derivatives at (i, j): x_der = (x_xi, y_xi),
    //                              y_der = (x_eta, y_eta)
    dgrid_ij.x_der.x = 1.0L;
    dgrid_ij.x_der.y = 2.0L;
    dgrid_ij.y_der.x = 3.0L;
    dgrid_ij.y_der.y = 4.0L;

    // Second derivatives at (i, j):
    //   x_der  = (x_xixi,   y_xixi)   -- unused for ip1_j stencil
    //   xy_der = (x_xieta,  y_xieta)
    //   y_der  = (x_etaeta, y_etaeta)
    d2grid_ij.x_der.x = 5.0L;
    d2grid_ij.x_der.y = 6.0L;
    d2grid_ij.xy_der.x = 7.0L;
    d2grid_ij.xy_der.y = 8.0L;
    d2grid_ij.y_der.x = 9.0L;
    d2grid_ij.y_der.y = 10.0L;


    // Closed-form expected values
    mult = 1.0L/dxi;
    expected_ddx_x = mult * ((coeff_ij.alpha/dxi)
                             + (dgrid_ij.x_der.x * d2grid_ij.y_der.x)
                             - (dgrid_ij.y_der.x * d2grid_ij.xy_der.x));
    expected_ddy_x = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x)
                             - (dgrid_ij.y_der.y * d2grid_ij.xy_der.x));
    expected_ddx_y = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y)
                             - (dgrid_ij.y_der.x * d2grid_ij.xy_der.y));
    expected_ddy_y = mult * ((coeff_ij.alpha/dxi)
                             + (dgrid_ij.x_der.y * d2grid_ij.y_der.y)
                             - (dgrid_ij.y_der.y * d2grid_ij.xy_der.y));


    // Call function
    ders = ell_interior_ip1_j (dxi, dgrid_ij, coeff_ij, d2grid_ij);


    // All four entries match closed form
    cr_assert_float_eq (ders.ddx.x, expected_ddx_x, 1e-15L,
                        "ddx.x mismatch: got %.15Le, expected %.15Le",
                        ders.ddx.x, expected_ddx_x);
    cr_assert_float_eq (ders.ddy.x, expected_ddy_x, 1e-15L,
                        "ddy.x mismatch: got %.15Le, expected %.15Le",
                        ders.ddy.x, expected_ddy_x);
    cr_assert_float_eq (ders.ddx.y, expected_ddx_y, 1e-15L,
                        "ddx.y mismatch: got %.15Le, expected %.15Le",
                        ders.ddx.y, expected_ddx_y);
    cr_assert_float_eq (ders.ddy.y, expected_ddy_y, 1e-15L,
                        "ddy.y mismatch: got %.15Le, expected %.15Le",
                        ders.ddy.y, expected_ddy_y);
}
