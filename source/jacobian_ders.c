#include <stdio.h>
#include "globvar.h"
#include "types.h"
#include "aux_functions.h"



/* Global variables */
long double           ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0, THREE = 3.0,
                      FOUR = 4.0, FIVE = 5.0, SIX = 6.0, SEVEN = 7.0, EIGHT = 8.0,
                      NINE = 9.0, TEN = 10.0;



/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for biharmonic grid generation   */
/*                                                           */
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------o------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*-----------------------------------------------------------*/
/*                Interior point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i - 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D interior_im1_jm1 (long double dxi, long double deta, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j - 1)
    ders.ddy.y                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j - 1)
    ders.ddp.p                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dp(x(i, j))/dp(i - 1, j - 1)
    ders.ddq.q                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dq(x(i, j))/dq(i - 1, j - 1)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i, j - 1)

   Input parameters: deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     pq_ij      - control function values (p, q) at point (i, j)
                     d2pq_ij    - second order derivatives of (p, q) wrt
                                  (xi, eta) at point (i, j)
*/
bh_jacobian_2D interior_i_jm1 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                               grid_dder_2D d2grid_ij, point_2D pq_ij, grid_dder_2D d2pq_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             cf_x, cf_y, temp, mult;


    /* Temporary variables */
    cf_x                                = (pq_ij.x * dgrid_ij.x_der.x) +
                                          (pq_ij.y * dgrid_ij.y_der.x);
    cf_y                                = (pq_ij.x * dgrid_ij.x_der.y) +
                                          (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.gamma/deta) - ((coeff_ij.J *
                                          coeff_ij.J * pq_ij.y)/TWO);
    mult                                = ONE/deta;

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) + (coeff_ij.J *
                                          dgrid_ij.x_der.y * cf_x) + temp); // dF(x(i, j))/dx(i, j - 1)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) + (coeff_ij.J *
                                          dgrid_ij.x_der.y * cf_y)); // dF(y(i, j))/dx(i, j - 1)
    ders.ddx.p                          = mult * (-(dgrid_ij.y_der.x * d2pq_ij.x_der.x) +
                                          (dgrid_ij.x_der.x * d2pq_ij.xy_der.x)); // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = mult * (-(dgrid_ij.y_der.x * d2pq_ij.x_der.y) +
                                          (dgrid_ij.x_der.x * d2pq_ij.xy_der.y)); // dF(q(i, j))/dx(i, j - 1)

    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) - (coeff_ij.J *
                                          dgrid_ij.x_der.x * cf_x)); // dF(x(i, j))/dy(i, j - 1)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) - (coeff_ij.J *
                                          dgrid_ij.x_der.x * cf_y) + temp); //dF(y(i, j))/dy(i, j - 1)
    ders.ddy.p                          = mult * (-(dgrid_ij.y_der.y * d2pq_ij.x_der.x) +
                                          (dgrid_ij.x_der.y * d2pq_ij.xy_der.x)); // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = mult * (-(dgrid_ij.y_der.y * d2pq_ij.x_der.y) +
                                          (dgrid_ij.x_der.y * d2pq_ij.xy_der.y)); // dF(q(i, j))/dy(i, j - 1)

    ders.ddp.p                          = coeff_ij.gamma/(deta * deta); // dF(p(i, j))/dp(i, j - 1)

    ders.ddq.q                          = coeff_ij.gamma/(deta * deta); // dF(q(i, j))/dq(i, j - 1)


    return ders;
}







/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i + 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D interior_ip1_jm1 (long double dxi, long double deta, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i + 1, j - 1)
    ders.ddy.y                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i + 1, j - 1)
    ders.ddp.p                          = (coeff_ij.beta)/(TWO * dxi * deta); // dp(x(i, j))/dp(i + 1, j - 1)
    ders.ddq.q                          = (coeff_ij.beta)/(TWO * dxi * deta); // dq(x(i, j))/dq(i + 1, j - 1)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i - 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     pq_ij      - control function values (p, q) at point (i, j)
                     d2pq_ij    - second order derivatives of (p, q) wrt
                                  (xi, eta) at point (i, j)
*/
bh_jacobian_2D interior_im1_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                   grid_dder_2D d2grid_ij, point_2D pq_ij, grid_dder_2D d2pq_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             cf_x, cf_y, temp, mult;


    /* Temporary variables */
    cf_x                                = (pq_ij.x * dgrid_ij.x_der.x) +
                                          (pq_ij.y * dgrid_ij.y_der.x);
    cf_y                                = (pq_ij.x * dgrid_ij.x_der.y) +
                                          (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.alpha/dxi) - ((coeff_ij.J *
                                          coeff_ij.J * pq_ij.x)/TWO);
    mult                                = ONE/dxi;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                          (dgrid_ij.x_der.x * d2grid_ij.y_der.x) - (coeff_ij.J *
                                          dgrid_ij.y_der.y * cf_x) + temp); // dF(x(i, j))/dx(i - 1, j)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                          (dgrid_ij.x_der.x * d2grid_ij.y_der.y) - (coeff_ij.J *
                                          dgrid_ij.y_der.y * cf_y)); // dF(y(i, j))/dx(i - 1, j)
    ders.ddx.p                          = mult * ((dgrid_ij.y_der.x * d2pq_ij.xy_der.x) -
                                          (dgrid_ij.x_der.x * d2pq_ij.y_der.x)); // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = mult * ((dgrid_ij.y_der.x * d2pq_ij.xy_der.y) -
                                          (dgrid_ij.x_der.x * d2pq_ij.y_der.y)); // dF(q(i, j))/dx(i - 1, j)

    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                          (dgrid_ij.x_der.y * d2grid_ij.y_der.x) + (coeff_ij.J *
                                          dgrid_ij.y_der.x * cf_x)); // dF(x(i, j))/dy(i - 1, j)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                          (dgrid_ij.x_der.y * d2grid_ij.y_der.y) + (coeff_ij.J *
                                          dgrid_ij.y_der.x * cf_y) + temp); // dF(y(i, j))/dy(i - 1, j)
    ders.ddy.p                          = mult * ((dgrid_ij.y_der.y * d2pq_ij.xy_der.x) -
                                          (dgrid_ij.x_der.y * d2pq_ij.y_der.x)); // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = mult * ((dgrid_ij.y_der.y * d2pq_ij.xy_der.y) -
                                          (dgrid_ij.x_der.y * d2pq_ij.y_der.y)); // dF(q(i, j))/dy(i - 1, j)

    ders.ddp.p                          = coeff_ij.alpha/(dxi * dxi); // dF(p(i, j))/dp(i - 1, j)

    ders.ddq.q                          = coeff_ij.alpha/(dxi * dxi); // dF(q(i, j))/dq(i - 1, j)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivative of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D interior_i_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(x(i, j))/dx(i, j)

    ders.ddy.y                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(y(i, j))/dy(i, j)

    ders.ddp.x                          = coeff_ij.J * coeff_ij.J * dgrid_ij.x_der.x;
    ders.ddp.y                          = coeff_ij.J * coeff_ij.J * dgrid_ij.x_der.y;
    ders.ddp.p                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(p(i, j))/dp(i, j)

    ders.ddq.x                          = coeff_ij.J * coeff_ij.J * dgrid_ij.y_der.x;
    ders.ddq.y                          = coeff_ij.J * coeff_ij.J * dgrid_ij.y_der.y;
    ders.ddq.q                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(p(i, j))/dp(i, j)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     pq_ij      - control function values (p, q) at point (i, j)
                     d2pq_ij    - second order derivatives of (p, q) wrt
                                  (xi, eta) at point (i, j)
*/
bh_jacobian_2D interior_ip1_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                   grid_dder_2D d2grid_ij, point_2D pq_ij, grid_dder_2D d2pq_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             cf_x, cf_y, temp, mult;


    /* Temporary variables */
    cf_x                                = (pq_ij.x * dgrid_ij.x_der.x) +
                                          (pq_ij.y * dgrid_ij.y_der.x);
    cf_y                                = (pq_ij.x * dgrid_ij.x_der.y) +
                                          (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.alpha/dxi) + ((coeff_ij.J *
                                          coeff_ij.J * pq_ij.x)/TWO);
    mult                                = ONE/dxi;

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                          (dgrid_ij.x_der.x * d2grid_ij.y_der.x) + (coeff_ij.J *
                                          dgrid_ij.y_der.y * cf_x) + temp); // dF(x(i, j))/dx(i + 1, j)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                          (dgrid_ij.x_der.x * d2grid_ij.y_der.y) + (coeff_ij.J *
                                          dgrid_ij.y_der.y * cf_y)); // dF(y(i, j))/dx(i + 1, j)
    ders.ddx.p                          = mult * (-(dgrid_ij.y_der.x * d2pq_ij.xy_der.x) +
                                          (dgrid_ij.x_der.x * d2pq_ij.y_der.x)); // dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = mult * (-(dgrid_ij.y_der.x * d2pq_ij.xy_der.y) +
                                          (dgrid_ij.x_der.x * d2pq_ij.y_der.y)); // dF(q(i, j))/dx(i + 1, j)

    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                          (dgrid_ij.x_der.y * d2grid_ij.y_der.x) - (coeff_ij.J *
                                          dgrid_ij.y_der.x * cf_x)); // dF(x(i, j))/dy(i + 1, j)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                          (dgrid_ij.x_der.y * d2grid_ij.y_der.y) - (coeff_ij.J *
                                          dgrid_ij.y_der.x * cf_y) + temp); // dF(y(i, j))/dy(i + 1, j)
    ders.ddy.p                          = mult * (-(dgrid_ij.y_der.y * d2pq_ij.xy_der.x) +
                                          (dgrid_ij.x_der.y * d2pq_ij.y_der.x)); // dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = mult * (-(dgrid_ij.y_der.y * d2pq_ij.xy_der.y) +
                                          (dgrid_ij.x_der.y * d2pq_ij.y_der.y)); // dF(q(i, j))/dy(i + 1, j)

    ders.ddp.p                          = coeff_ij.alpha/(dxi * dxi); // dF(p(i, j))/dp(i + 1, j)

    ders.ddq.q                          = coeff_ij.alpha/(dxi * dxi); // dF(q(i, j))/dq(i + 1, j)


    return ders;
}







/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i - 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D interior_im1_jp1 (long double dxi, long double deta, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j + 1)
    ders.ddy.y                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j + 1)
    ders.ddp.p                          = (coeff_ij.beta)/(TWO * dxi * deta); // dp(x(i, j))/dp(i - 1, j + 1)
    ders.ddq.q                          = (coeff_ij.beta)/(TWO * dxi * deta); // dq(x(i, j))/dq(i - 1, j + 1)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i, j + 1)

   Input parameters: deta        - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     pq_ij      - control function values (p, q) at point (i, j)
                     d2pq_ij    - second order derivatives of (p, q) wrt
                                  (xi, eta) at point (i, j)
*/
bh_jacobian_2D interior_i_jp1 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                   grid_dder_2D d2grid_ij, point_2D pq_ij, grid_dder_2D d2pq_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             cf_x, cf_y, temp, mult;


    /* Temporary variables */
    cf_x                                = (pq_ij.x * dgrid_ij.x_der.x) +
                                          (pq_ij.y * dgrid_ij.y_der.x);
    cf_y                                = (pq_ij.x * dgrid_ij.x_der.y) +
                                          (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.gamma/deta) + ((coeff_ij.J *
                                          coeff_ij.J * pq_ij.x)/TWO);
    mult                                = ONE/deta;

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) - (coeff_ij.J *
                                          dgrid_ij.x_der.y * cf_x) + temp); // dF(x(i, j))/dx(i, j + 1)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) - (coeff_ij.J *
                                          dgrid_ij.x_der.y * cf_y)); // dF(y(i, j))/dx(i, j + 1)
    ders.ddx.p                          = mult * ((dgrid_ij.y_der.x * d2pq_ij.x_der.x) -
                                          (dgrid_ij.x_der.x * d2pq_ij.xy_der.x)); // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = mult * ((dgrid_ij.y_der.x * d2pq_ij.x_der.y) -
                                          (dgrid_ij.x_der.x * d2pq_ij.xy_der.y)); // dF(q(i, j))/dx(i, j + 1)

    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) + (coeff_ij.J *
                                          dgrid_ij.x_der.x * cf_x)); // dF(x(i, j))/dy(i, j + 1)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) + (coeff_ij.J *
                                          dgrid_ij.x_der.x * cf_y) + temp); //dF(y(i, j))/dy(i, j + 1)
    ders.ddy.p                          = mult * ((dgrid_ij.y_der.y * d2pq_ij.x_der.x) -
                                          (dgrid_ij.x_der.y * d2pq_ij.xy_der.x)); // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = mult * ((dgrid_ij.y_der.y * d2pq_ij.x_der.y) -
                                          (dgrid_ij.x_der.y * d2pq_ij.xy_der.y)); // dF(q(i, j))/dy(i, j + 1)

    ders.ddp.p                          = coeff_ij.gamma/(deta * deta); // dF(p(i, j))/dp(i, j + 1)

    ders.ddq.q                          = coeff_ij.gamma/(deta * deta); // dF(q(i, j))/dq(i, j + 1)


    return ders;
}






/*
   Compute derivatives of interior F(u(i, j))
   with respect to u(i + 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D interior_ip1_jp1 (long double dxi, long double deta, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i + 1, j + 1)
    ders.ddy.y                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i + 1, j + 1)
    ders.ddp.p                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dp(x(i, j))/dp(i + 1, j + 1)
    ders.ddq.q                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dq(x(i, j))/dq(i + 1, j + 1)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for interior
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     pq_ij      - values of control functions (p, q) at
                                  point (i1, j1)
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     d2pq_ij    - values of second-order derivatives of
                                  (p, q) at point (i1, j1)
*/
bh_jacobian_2D interior_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                    long double deta, point_2D pq_ij, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                                    grid_dder_2D d2pq_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 - 1 && j2 == j1 - 1)   // dF(u(i1, j1))/du(i1 - 1, j1 - 1)
    {
        dFdu                            = interior_im1_jm1 (dxi, deta, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = interior_i_jm1 (deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, pq_ij, d2pq_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 1, j1 - 1)
    {
        dFdu                            = interior_ip1_jm1 (dxi, deta, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = interior_im1_j (dxi, dgrid_ij, coeff_ij,
                                                          d2grid_ij, pq_ij, d2pq_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = interior_i_j (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = interior_ip1_j (dxi, dgrid_ij, coeff_ij,
                                                          d2grid_ij, pq_ij, d2pq_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 1, j1 + 1)
    {
        dFdu                            = interior_im1_jp1 (dxi, deta, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = interior_i_jp1 (deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, pq_ij, d2pq_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 1, j1 + 1)
    {
        dFdu                            = interior_ip1_jp1 (dxi, deta, coeff_ij);
    }
    else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      o------x------x------x               */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*-----------------------------------------------------------*/
/*                  xi = 0 point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta        - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi0_i_jm1 (long double dxi, long double deta,
                          grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                          point_2D Dx_Dy_ij, point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -(coeff_ij.gamma/deta) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));


    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 1)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_ip1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                            coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 2, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_ip2_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                            coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 2, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 2, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_i_j (long double dxi, long double deta, long double c,
                        grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                        point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variable */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = (TWO/THREE) * (-(coeff_ij.alpha/dxi) +
                                          ((coeff_ij.gamma * dxi)/(deta * deta)));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + ((ONEHALF * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + ((ONEHALF * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (ONEHALF/dxi) * ((nablas_ij.y * dgrid_ij.y_der.y) -
                                                           Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (ONEHALF/dxi) * ((nablas_ij.y * dgrid_ij.y_der.x) -
                                                           Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE; //dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx; // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx + ((ONEHALF * c * n_ij.y)/dxi); // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE; // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy; // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady - ((ONEHALF * c * n_ij.x)/dxi); // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE; // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE; // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_ip1_j (long double dxi, long double c, grid_der_2D dgrid_ij,
                          coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                          point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = (FIVE * coeff_ij.alpha)/(FOUR * dxi);

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                       Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                       Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx; // dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = -detadx - ((TWO * c * n_ij.y)/dxi); // dF(q(i, j))/dx(i + 1, j)
    ders.ddy.p                          = -dxidy; // dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = -detady + ((TWO * c * n_ij.x)/dxi); // dF(q(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_ip2_j (long double dxi, long double c, grid_der_2D dgrid_ij,
                          coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                          point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(FOUR * coeff_ij.alpha)/dxi;

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx; // dF(p(i, j))/dx(i + 2, j)
    ders.ddx.q                          = -detadx + ((HALF * c * n_ij.y)/dxi); // dF(q(i, j))/dx(i + 2, j)
    ders.ddy.p                          = -dxidy; // dF(p(i, j))/dy(i + 2, j)
    ders.ddy.q                          = -detady - ((HALF * c * n_ij.x)/dxi); // dF(q(i, j))/dy(i + 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 3, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_ip3_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.alpha)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * dxi);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 3, j)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi0_i_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                              coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                              point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;

    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -(coeff_ij.gamma/deta) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));


    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_ip1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 0 F(u(i, j))
   with respect to u(i + 2, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_ip2_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 2, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 2, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j + 1)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 0
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi0_point_ders (int i1, int j1, int i2, int j2,
                               long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                               grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                               point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 && j2 == j1 - 1)   // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = xi0_i_jm1 (dxi, deta, dgrid_ij, coeff_ij,
                                                     d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 1, j1 - 1)
    {
        dFdu                            = xi0_ip1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 2, j1 - 1)
    {
        dFdu                            = xi0_ip2_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi0_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                   d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = xi0_ip1_j (dxi, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                     Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 2, j1)
    {
        dFdu                            = xi0_ip2_j (dxi, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                     Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 3 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 3, j1)
    {
        dFdu                            = xi0_ip3_j (dxi, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = xi0_i_jp1 (dxi, deta, dgrid_ij, coeff_ij,
                                                     d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 1, j1 + 1)
    {
        dFdu                            = xi0_ip1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 2, j1 + 1)
    {
        dFdu                            = xi0_ip2_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                      x------x------x                      */
/*                      |      |      |                      */
/*               x------x------x------o                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*-----------------------------------------------------------*/
/*                  xi = 1 point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 2, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_im2_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 2, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 2, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_im1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi1_i_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                              coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                              point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -(coeff_ij.gamma/deta) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + (nablas_ij.x/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + (nablas_ij.x/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 1)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 3, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_im3_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.alpha)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * dxi);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 3, j)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_im2_j (long double dxi, long double c, grid_der_2D dgrid_ij,
                          coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                          point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(FOUR * coeff_ij.alpha)/dxi;

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 2, j)
    ders.ddx.q                          = -detadx - ((HALF * c * n_ij.y)/dxi);  // dF(q(i, j))/dx(i - 2, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 2, j)
    ders.ddy.q                          = -detady + ((HALF * c * n_ij.x)/dxi);  // dF(q(i, j))/dy(i - 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_im1_j (long double dxi, long double c, grid_der_2D dgrid_ij,
                          coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                          point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = (FIVE * coeff_ij.alpha)/(FOUR * dxi);

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) + temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) + temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                       Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                       Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = -detadx + ((TWO * c *n_ij.y)/dxi);  // dF(q(i, j))/dx(i - 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = -detady - ((TWO * c * n_ij.x)/dxi);  // dF(q(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_i_j (long double dxi, long double deta, long double c, grid_der_2D dgrid_ij,
                        coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                        point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variable */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(TWO/THREE) * ((coeff_ij.alpha/dxi) -
                                          ((coeff_ij.gamma * dxi)/(deta * deta)));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) + temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) + temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - ((ONEHALF * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - ((ONEHALF * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (ONEHALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                           Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (ONEHALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                           Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx - ((ONEHALF * c * n_ij.y)/dxi);  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady + ((ONEHALF * c * n_ij.x)/dxi);  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 2, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_im2_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 2, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 2, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i - 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_im1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 1 F(u(i, j))
   with respect to u(i, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi1_i_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                              coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                              point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J  * coeff_ij.J * deta);
    temp                                = -(coeff_ij.gamma/deta) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi1_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                               long double deta, long double c, grid_der_2D dgrid_ij,
                               coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                               point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 - 2 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 - 2, j1 - 1)
    {
        dFdu                            = xi1_im2_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 - 1, j1 - 1)
    {
        dFdu                            = xi1_im1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 1)   // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = xi1_i_jm1 (dxi, deta, dgrid_ij, coeff_ij,
                                                     d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 - 3 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 3, j1)
    {
        dFdu                            = xi1_im3_j (dxi, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 2, j1)
    {
        dFdu                            = xi1_im2_j (dxi, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                     Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = xi1_im1_j (dxi, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                     Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi1_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                   d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 2 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 2, j1 + 1)
    {
        dFdu                            = xi1_im2_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 1, j1 + 1)
    {
        dFdu                            = xi1_im1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = xi1_i_jp1 (dxi, deta, dgrid_ij, coeff_ij,
                                                     d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                             x                             */
/*                             |                             */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------o------x                      */
/*-----------------------------------------------------------*/
/*                 eta = 0 point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i - 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D eta0_im1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(coeff_ij.alpha/dxi) +
                                           ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) + temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) + temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                         Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                         Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D eta0_i_j (long double dxi, long double deta, long double c,
                         grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                         point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variable */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = (TWO/THREE) * (-(coeff_ij.gamma/deta) +
                                          ((coeff_ij.alpha * deta)/(dxi * dxi)));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) + temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) + temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (ONEHALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                            Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (ONEHALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                            Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + ((ONEHALF * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + ((ONEHALF * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx - ((ONEHALF * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy + ((ONEHALF * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D eta0_ip1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(coeff_ij.alpha/dxi) -
                                           ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) + temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) + temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i + 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i - 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta0_im1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j + 1)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta0_i_jp1 (long double deta, long double c, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = (FIVE * coeff_ij.gamma)/(FOUR * deta);

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                        Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                        Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy - ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i + 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta0_ip1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j + 1)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i - 1, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta0_im1_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j + 2)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta0_i_jp2 (long double deta, long double c, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -(FOUR * coeff_ij.gamma)/deta;

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


   /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j + 2)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 2)
    ders.ddy.p                          = -dxidy + ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j + 2)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 2)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i + 1, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta0_ip1_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j + 2)


    return ders;
}






/*
   Compute derivatives of eta = 0 F(u(i, j))
   with respect to u(i, j + 3)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta0_i_jp3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.gamma)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j + 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i, j + 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i, j + 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j + 3)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D eta0_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                long double deta, long double c, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 - 1 && j2 == j1)   // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = eta0_im1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                      d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = eta0_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                    d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = eta0_ip1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                      d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 1, j1 + 1)
    {
        dFdu                            = eta0_im1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = eta0_i_jp1 (deta, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                      Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 1, j1 + 1)
    {
        dFdu                            = eta0_ip1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 - 1, j1 + 2)
    {
        dFdu                            = eta0_im1_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1, j1 + 2)
    {
        dFdu                            = eta0_i_jp2 (deta, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                      Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 + 1, j1 + 2)
    {
        dFdu                            = eta0_ip1_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 3)    // dF(u(i1, j1))/du(i1, j1 + 3)
    {
        dFdu                            = eta0_i_jp3 (deta, dgrid_ij, coeff_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}








/*-----------------------------------------------------------*/
/*                      x------o------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                             |                             */
/*                             x                             */
/*-----------------------------------------------------------*/
/*                  eta = 1 point derivatives                */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i, j - 3)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta1_i_jm3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.gamma)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j + 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i, j + 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i, j + 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j + 3)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i - 1, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta1_im1_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j - 2)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta1_i_jm2 (long double deta, long double c, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -(FOUR * coeff_ij.gamma)/deta;

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j - 2)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 2)
    ders.ddy.p                          = -dxidy - ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j - 2)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 2)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i + 1, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta1_ip1_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j - 2)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i - 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta1_im1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i - 1, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i - 1, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta1_i_jm1 (long double deta, long double c, grid_der_2D dgrid_ij,
                           coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                           point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = (FIVE * coeff_ij.gamma)/(FOUR * deta);

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                        Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                        Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 1)
    ders.ddy.p                          = -dxidy + ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i + 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D eta1_ip1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                 coeffs_1 coeff_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;
    long double             temp;


    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dy(i + 1, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dx(i + 1, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j - 1)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta1_im1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                               coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                               point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(coeff_ij.alpha/dxi) -
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta1_i_j (long double dxi, long double deta, long double c, grid_der_2D dgrid_ij,
                         coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                         point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variable */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = (TWO/THREE) * (((coeff_ij.alpha * deta)/(dxi *
                                          dxi)) - (coeff_ij.gamma/deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (ONEHALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                            Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (ONEHALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                            Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - ((ONEHALF * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - ((ONEHALF * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx + ((ONEHALF * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy - ((ONEHALF * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of eta = 1 F(u(i, j))
   with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D eta1_ip1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                               coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                               point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -(coeff_ij.alpha/dxi) +
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i + 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D eta1_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                long double deta, long double c, grid_der_2D dgrid_ij,
                                coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 && j2 == j1 - 3)   // dF(u(i1, j1))/du(i1, j1 - 3)
    {
        dFdu                            = eta1_i_jm3 (deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 - 1, j1 - 2)
    {
        dFdu                            = eta1_im1_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1, j1 - 2)
    {
        dFdu                            = eta1_i_jm2 (deta, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                      Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 + 1, j1 - 2)
    {
        dFdu                            = eta1_ip1_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 - 1, j1 - 1)
    {
        dFdu                            = eta1_im1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = eta1_i_jm1 (deta, c, dgrid_ij, coeff_ij, d2grid_ij,
                                                      Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 1, j1 - 1)
    {
        dFdu                            = eta1_ip1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = eta1_im1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                      d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = eta1_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                    d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = eta1_ip1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                      d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }

    return dFdu;
}








/*-----------------------------------------------------------*/
/*                      x                                    */
/*                      |                                    */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      o------x------x------x               */
/*-----------------------------------------------------------*/
/*             xi = 0, eta = 0 point derivatives             */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta0_i_j (long double dxi, long double deta, long double c,
                             grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                             point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult, comp1, comp2;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variables */
    temp                               = (TWO/THREE) * (((coeff_ij.alpha * deta)/dxi) +
                                          ((coeff_ij.gamma * dxi)/deta));
    mult                               = THREE/(coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Derivatives of Dx and Dy */
    comp1                               = (dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x);
    dDxdx                               = mult * ((comp1 * dxi) + (comp2 * deta) - temp +
                                                  (ONEHALF * coeff_ij.beta));
    comp1                               = (dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x);
    dDxdy                               = mult * ((comp1 * dxi) + (comp2 * deta));
    comp1                               = (dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.x* d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y);
    dDydx                               = mult * ((comp1 * dxi) + (comp2 * deta));
    comp1                               = (dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y);
    dDydy                               = mult * ((comp1 * dxi) + (comp2 * deta) - temp +
                                                  (ONEHALF * coeff_ij.beta));

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.y/dxi) -
                                          (dgrid_ij.x_der.y/deta))) + (Dx_Dy_ij.y/deta);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.x * ((dgrid_ij.x_der.x/deta) -
                                          (dgrid_ij.y_der.x/dxi))) - (Dx_Dy_ij.x/deta);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.y/dxi) -
                                          (dgrid_ij.x_der.y/deta))) - (Dx_Dy_ij.y/dxi);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.x_der.x/deta) -
                                          (dgrid_ij.y_der.x/dxi))) + (Dx_Dy_ij.x/dxi);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) +
                                                  (ONEHALF * temp));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx - ((ONEHALF * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy + ((ONEHALF * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = ((FIVE * coeff_ij.alpha)/(FOUR * dxi)) -
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                          Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                          Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i + 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip2_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variable */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -((FOUR * coeff_ij.alpha)/dxi) +
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                          Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                          Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i + 2, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i + 2, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i + 2, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i + 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 3, j)

   Input parameters: deta       - spacing in eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip3_j (long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.alpha)/(coeff_ij.J * coeff_ij.J *
                                          coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 3, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta0_i_jp1 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = ((FIVE * coeff_ij.gamma)/(FOUR * deta)) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                          Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                          Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy - ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(EIGHT * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 1, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 1, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 2, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip2_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 2, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 2, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta0_i_jp2 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -((FOUR * coeff_ij.gamma)/deta) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j + 2)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 2)
    ders.ddy.p                          = -dxidy + ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j + 2)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 1, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip1_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 1, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 1, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j + 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i + 2, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_ip2_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 2, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 2, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j + 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 0
   F(u(i, j)) with respect to u(i, j + 3)

   Input parameters: deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta0_i_jp3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.gamma)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j + 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i, j + 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i, j + 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j + 3)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi0_eta0_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                    long double deta, long double c, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                    point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 && j2 == j1)   // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi0_eta0_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                        d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = xi0_eta0_ip1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 2, j1)
    {
        dFdu                            = xi0_eta0_ip2_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 3 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 3, j1)
    {
        dFdu                            = xi0_eta0_ip3_j (dxi, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = xi0_eta0_i_jp1 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 1, j1 + 1)
    {
        dFdu                            = xi0_eta0_ip1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 + 2, j1 + 1)
    {
        dFdu                            = xi0_eta0_ip2_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1, j1 + 2)
    {
        dFdu                            = xi0_eta0_i_jp2 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 + 1, j1 + 2)
    {
        dFdu                            = xi0_eta0_ip1_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 + 2, j1 + 2)
    {
        dFdu                            = xi0_eta0_ip2_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 3)    // dF(u(i1, j1))/du(i1, j1 + 3)
    {
        dFdu                            = xi0_eta0_i_jp3 (deta, dgrid_ij, coeff_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                                    x                      */
/*                                    |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*               x------x------x------o                      */
/*-----------------------------------------------------------*/
/*             xi = 1, eta = 0 point derivatives             */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 3, j)

   Input parameters: deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im3_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.alpha)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * dxi);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 3, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im2_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -((FOUR * coeff_ij.alpha)/dxi) -
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 2, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 2, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 2, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = ((FIVE * coeff_ij.alpha)/(FOUR * dxi)) +
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta0_i_j (long double dxi, long double deta, long double c,
                             grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                             point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult, comp1, comp2;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variables */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * dxi * deta);
    temp                                = (TWO/THREE) * (((coeff_ij.alpha * deta)/dxi) +
                                                         ((coeff_ij.gamma * dxi)/deta));

    /* Derivatives of Dx and Dy */
    comp1                               = (dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x);
    dDxdx                               = mult * ((comp1 * dxi) - (comp2 * deta) -
                                                  temp - (ONEHALF * coeff_ij.beta));

    comp1                               = (dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x);
    dDxdy                               = mult * ((comp1 * dxi) - (comp2 * deta));

    comp1                               = (dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y);
    dDydx                               = mult * ((comp1 * dxi) - (comp2 * deta));

    comp1                               = (dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y);
    dDydy                               = mult * ((comp1 * dxi) - (comp2 * deta) -
                                                  temp - (ONEHALF * coeff_ij.beta));

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.y/dxi) +
                                          (dgrid_ij.x_der.y/deta))) - (Dx_Dy_ij.y/deta);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.x/dxi) +
                                          (dgrid_ij.x_der.x/deta))) - (Dx_Dy_ij.x/deta);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.y/dxi) +
                                          (dgrid_ij.x_der.y/deta))) - (Dx_Dy_ij.y/dxi);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.x/dxi) +
                                          (dgrid_ij.x_der.x/deta))) - (Dx_Dy_ij.x/dxi);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) +
                                                  (ONEHALF * temp));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx - ((ONEHALF * c * n_ij.y)/deta); // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(p(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy + ((ONEHALF * c * n_ij.x)/deta); // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 2, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im2_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 2, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 2, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 1, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im1_jp1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (EIGHT * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j + 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 2, j + 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 2, j + 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j + 1)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i, j + 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta0_i_jp1 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = ((FIVE * coeff_ij.gamma)/(FOUR * deta)) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                        Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                        Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy - ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)

    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 2, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im2_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 2, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 2, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j + 2)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i - 1, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_im1_jp2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j + 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 1, j + 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 1, j + 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j + 2)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i, j + 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta0_i_jp2 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -((FOUR * coeff_ij.gamma)/deta) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j + 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j + 1)
    ders.ddy.p                          = -dxidy + ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j + 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j + 1)

    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 0
   F(u(i, j)) with respect to u(i, j + 3)

   Input parameters: deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta0_i_jp3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.gamma)/(coeff_ij.J * coeff_ij.J *
                                          coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j + 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i, j + 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i, j + 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j + 3)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi1_eta0_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                    long double deta, long double c, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                    point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 - 3 && j2 == j1)   // dF(u(i1, j1))/du(i1 - 3, j1)
    {
        dFdu                            = xi1_eta0_im3_j (dxi, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 2, j1)
    {
        dFdu                            = xi1_eta0_im2_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 - 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = xi1_eta0_im1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi1_eta0_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                        d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 2 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 2, j1 + 1)
    {
        dFdu                            = xi1_eta0_im2_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1 - 1, j1 + 1)
    {
        dFdu                            = xi1_eta0_im1_jp1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 1)    // dF(u(i1, j1))/du(i1, j1 + 1)
    {
        dFdu                            = xi1_eta0_i_jp1 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 2 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 - 2, j1 + 2)
    {
        dFdu                            = xi1_eta0_im2_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1 - 1, j1 + 2)
    {
        dFdu                            = xi1_eta0_im1_jp2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 + 2)    // dF(u(i1, j1))/du(i1, j1 + 2)
    {
        dFdu                            = xi1_eta0_i_jp2 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 && j2 == j1 + 3)    // dF(u(i1, j1))/du(i1, j1 + 3)
    {
        dFdu                            = xi1_eta0_i_jp3 (deta, dgrid_ij, coeff_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                      o------x------x------x               */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |                                    */
/*                      x                                    */
/*-----------------------------------------------------------*/
/*             xi = 0, eta = 1 point derivatives             */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i, j - 3)

   Input parameters: deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_i_jm3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.gamma)/(coeff_ij.J * coeff_ij.J *
                                          coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j - 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i, j - 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i, j - 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j - 3)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta1_i_jm2 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -((FOUR * coeff_ij.gamma)/deta) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j - 2)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 2)
    ders.ddy.p                          = -dxidy - ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j - 2)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 1, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip1_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 1, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 1, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 1, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 1, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 2, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip2_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 2, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 2, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta1_i_jm1 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = ((FIVE * coeff_ij.gamma)/(FOUR * deta)) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                        Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                        Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 1)
    ders.ddy.p                          = -dxidy + ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (EIGHT * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 2, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 2, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 2, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip2_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 2, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 2, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 2, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 2, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta1_i_j (long double dxi, long double deta, long double c,
                             grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                             point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult, comp1, comp2;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variables */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * dxi * deta);
    temp                                = (TWO/THREE) * (((coeff_ij.alpha * deta)/dxi) +
                                                         ((coeff_ij.gamma * dxi)/deta));

    /* Derivatives of Dx and Dy */
    comp1                               = -(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x);
    dDxdx                               = mult * ((comp1 * dxi) + (comp2 * deta) -
                                                  temp - (ONEHALF * coeff_ij.beta));

    comp1                               = -(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x);
    dDxdy                               = mult * ((comp1 * dxi) + (comp2 * deta));

    comp1                               = -(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y);
    dDydx                               = mult * ((comp1 * dxi) + (comp2 * deta));

    comp1                               = -(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y);
    dDydy                               = mult * ((comp1 * dxi) + (comp2 * deta) -
                                                  temp - (ONEHALF * coeff_ij.beta));

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.y/dxi) +
                                          (dgrid_ij.x_der.y/deta))) - (Dx_Dy_ij.y/deta);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.x/dxi) +
                                          (dgrid_ij.x_der.x/deta))) - (Dx_Dy_ij.x/deta);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.y/dxi) +
                                          (dgrid_ij.x_der.y/deta))) - (Dx_Dy_ij.y/dxi);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) +
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.x/dxi) +
                                          (dgrid_ij.x_der.x/deta))) - (Dx_Dy_ij.x/dxi);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) -
                                                  (ONEHALF * temp));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx + ((ONEHALF * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy - ((ONEHALF * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = ((FIVE * coeff_ij.alpha)/(FOUR * dxi)) +
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                       Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                       Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   //  dF(p(i, j))/dx(i + 1, j)
    ders.ddx.q                          = -detadx;  //  dF(q(i, j))/dx(i + 1, j)
    ders.ddy.p                          = -dxidy;   //  dF(p(i, j))/dy(i + 1, j)
    ders.ddy.q                          = -detady;  //  dF(q(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip2_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -((FOUR * coeff_ij.alpha)/dxi) -
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   //  dF(p(i, j))/dx(i + 2, j)
    ders.ddx.q                          = -detadx;  //  dF(q(i, j))/dx(i + 2, j)
    ders.ddy.p                          = -dxidy;   //  dF(p(i, j))/dy(i + 2, j)
    ders.ddy.q                          = -detady;  //  dF(q(i, j))/dy(i + 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 0, eta = 1
   F(u(i, j)) with respect to u(i + 3, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi0_eta1_ip3_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.alpha)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * dxi);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i + 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i + 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i + 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i + 3, j)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi0_eta1_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                    long double deta, long double c, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                    point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 && j2 == j1 - 3)   // dF(u(i1, j1))/du(i1, j1 - 3)
    {
        dFdu                            = xi0_eta1_i_jm3 (deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1, j1 - 2)
    {
        dFdu                            = xi0_eta1_i_jm2 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 + 1, j1 - 2)
    {
        dFdu                            = xi0_eta1_ip1_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 + 2, j1 - 2)
    {
        dFdu                            = xi0_eta1_ip2_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = xi0_eta1_i_jm1 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 1, j1 - 1)
    {
        dFdu                            = xi0_eta1_ip1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 + 2 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 + 2, j1 - 1)
    {
        dFdu                            = xi0_eta1_ip2_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi0_eta1_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                        d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 + 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 1, j1)
    {
        dFdu                            = xi0_eta1_ip1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 2, j1)
    {
        dFdu                            = xi0_eta1_ip2_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 + 3 && j2 == j1)    // dF(u(i1, j1))/du(i1 + 3, j1)
    {
        dFdu                            = xi0_eta1_ip3_j (dxi, dgrid_ij, coeff_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*               x------x------x------o                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*                                    |                      */
/*                                    x                      */
/*-----------------------------------------------------------*/
/*             xi = 1, eta = 1 point derivatives             */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i, j - 3)

   Input parameters: deta       - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_i_jm3 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.gamma)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * deta * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i, j - 3)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i, j - 3)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i, j - 3)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i, j - 3)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 2, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im2_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.beta)/(TWO * coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 2, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 2, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 1, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im1_jm2 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j - 2)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 1, j - 2)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 1, j - 2)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i, j - 2)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta1_i_jm2 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = -((FOUR * coeff_ij.gamma)/deta) +
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (HALF/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                         Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) + temp);
    temp                                = (HALF/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                         Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) - temp);
    temp                                = dDxdx - (nablas_ij.y/(TWO * deta));
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy - (nablas_ij.y/(TWO * deta));
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx + ((c * n_ij.y)/(TWO * deta));   // dF(p(i, j))/dx(i, j - 2)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 2)
    ders.ddy.p                          = -dxidy - ((c * n_ij.x)/(TWO * deta));   // dF(p(i, j))/dy(i, j - 2)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 2)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 2, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im2_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = (TWO * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 2, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 2, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 2, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 2, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 1, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im1_jm1 (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                     coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(EIGHT * coeff_ij.beta)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * deta);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 1, j - 1)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 1, j - 1)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 1, j - 1)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i, j - 1)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta1_i_jm1 (long double dxi, long double deta, long double c,
                               grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                               point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * deta);
    temp                                = ((FIVE * coeff_ij.gamma)/(FOUR * deta)) -
                                          ((THREE * coeff_ij.beta)/(TWO * dxi));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (TWO/deta) * ((dgrid_ij.x_der.y * nablas_ij.x) -
                                                        Dx_Dy_ij.y);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) - temp);
    temp                                = (TWO/deta) * ((dgrid_ij.x_der.x * nablas_ij.x) -
                                                        Dx_Dy_ij.x);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) + temp);
    temp                                = dDxdx + ((TWO * nablas_ij.y)/deta);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * temp));
    temp                                = dDydy + ((TWO * nablas_ij.y)/deta);
    detady                              = mult * ((dgrid_ij.x_der.x * temp) -
                                                  (dgrid_ij.x_der.y * dDxdy));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx - ((TWO * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j - 1)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j - 1)
    ders.ddy.p                          = -dxidy + ((TWO * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j - 1)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 3, j)

   Input parameters: deta       - spacing of xi in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im3_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij)
{
    bh_jacobian_2D          ders;
    long double             temp;


    /* Temporary variable */
    temp                                = -(coeff_ij.alpha)/(coeff_ij.J *
                                          coeff_ij.J * coeff_ij.J * dxi * dxi);

    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = temp * dgrid_ij.y_der.y;  // dF(p(i, j))/dx(i - 3, j)
    ders.ddx.q                          = -temp * dgrid_ij.x_der.y; // dF(q(i, j))/dx(i - 3, j)
    ders.ddy.p                          = -temp * dgrid_ij.y_der.x; // dF(p(i, j))/dy(i - 3, j)
    ders.ddy.q                          = temp * dgrid_ij.x_der.x;  // dF(q(i, j))/dy(i - 3, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 2, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im2_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = ONE/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = -((FOUR * coeff_ij.alpha)/dxi) +
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * (-(dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * (-(dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx - (nablas_ij.x/(TWO * dxi));
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy - (nablas_ij.x/(TWO * dxi));
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                        Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) - temp);
    temp                                = (HALF/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) + temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 2, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 2, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 2, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 2, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i - 1, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta1_im1_j (long double dxi, long double deta, grid_der_2D dgrid_ij,
                                   coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                   point_2D nablas_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult;


    /* Temporary variables */
    mult                                = FOUR/(coeff_ij.J * coeff_ij.J * dxi);
    temp                                = ((FIVE * coeff_ij.alpha)/(FOUR * dxi)) -
                                          ((THREE * coeff_ij.beta)/(TWO * deta));

    /* Derivatives of Dx and Dy */
    dDxdx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x) +
                                                  temp);
    dDxdy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x));
    dDydx                               = mult * ((dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y));
    dDydy                               = mult * ((dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y) +
                                                  temp);

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = dDxdx + ((TWO * nablas_ij.x)/dxi);
    dxidx                               = mult * ((dgrid_ij.y_der.y * temp) -
                                                  (dgrid_ij.y_der.x * dDydx));
    temp                                = dDydy + ((TWO * nablas_ij.x)/dxi);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * temp));
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.y * nablas_ij.y) -
                                                       Dx_Dy_ij.y);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) + temp);
    temp                                = (TWO/dxi) * ((dgrid_ij.y_der.x * nablas_ij.y) -
                                                        Dx_Dy_ij.x);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) - temp);


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.p                          = -dxidx;   // dF(p(i, j))/dx(i - 1, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i - 1, j)
    ders.ddy.p                          = -dxidy;   // dF(p(i, j))/dy(i - 1, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of xi = 1, eta = 1
   F(u(i, j)) with respect to u(i, j)

   Input parameters: dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - first derivatives of (x, y) wrt (xi, eta)
                                  at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) wrt
                                  (xi, eta) at point (i, j)
                     Dx_Dy_ij   - Dx and Dy at point (i, j)
                     nablas_ij  - Laplacians of xi and eta at point (i, j)
                     n_ij       - outward normal vector at point (i, j)
*/
bh_jacobian_2D xi1_eta1_i_j (long double dxi, long double deta, long double c,
                             grid_der_2D dgrid_ij, coeffs_1 coeff_ij, grid_dder_2D d2grid_ij,
                             point_2D Dx_Dy_ij, point_2D nablas_ij, point_2D n_ij)
{
    /* Return ders */
    bh_jacobian_2D          ders;

    /* Local variables */
    long double             dDxdx, dDxdy, dDydx, dDydy,
                            dxidx, dxidy, detadx, detady,
                            temp, mult, comp1, comp2;
    const long double       ONEHALF = THREE/TWO;


    /* Temporary variables */
    mult                                = THREE/(coeff_ij.J * coeff_ij.J * dxi * deta);
    temp                                = (TWO/THREE) * (((coeff_ij.alpha * deta)/dxi) +
                                                         ((coeff_ij.gamma * dxi)/deta));

    /* Derivatives of Dx and Dy */
    comp1                               = -(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.x);
    dDxdx                               = mult * ((comp1 * dxi) - (comp2 * deta) -
                                                  temp + (ONEHALF * coeff_ij.beta));

    comp1                               = -(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.x);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.x);
    dDxdy                               = mult * ((comp1 * dxi) - (comp2 * deta));

    comp1                               = -(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.x_der.y * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.y_der.y * Dx_Dy_ij.y);
    dDydx                               = mult * ((comp1 * dxi) - (comp2 * deta));

    comp1                               = -(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                          (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                          (coeff_ij.J * dgrid_ij.x_der.x * Dx_Dy_ij.y);
    comp2                               = (dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                          (dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                          (coeff_ij.J * dgrid_ij.y_der.x * Dx_Dy_ij.y);
    dDydy                               = mult * ((comp1 * dxi) - (comp2 * deta) -
                                                  temp + (ONEHALF * coeff_ij.beta));

    /* Derivatives of Laplacians of xi and eta */
    mult                                = ONE/coeff_ij.J;

    temp                                = (nablas_ij.x * ((dgrid_ij.y_der.y)/dxi) -
                                          (dgrid_ij.x_der.y/deta)) + (Dx_Dy_ij.y/deta);
    dxidx                               = mult * ((dgrid_ij.y_der.y * dDxdx) -
                                                  (dgrid_ij.y_der.x * dDydx) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.x * ((dgrid_ij.x_der.x)/deta) -
                                          (dgrid_ij.y_der.x/dxi)) - (Dx_Dy_ij.x/deta);
    dxidy                               = mult * ((dgrid_ij.y_der.y * dDxdy) -
                                                  (dgrid_ij.y_der.x * dDydy) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.y_der.y)/dxi) -
                                          (dgrid_ij.x_der.y/deta)) - (Dx_Dy_ij.y/dxi);
    detadx                              = mult * ((dgrid_ij.x_der.x * dDydx) -
                                                  (dgrid_ij.x_der.y * dDxdx) -
                                                  (ONEHALF * temp));
    temp                                = (nablas_ij.y * ((dgrid_ij.x_der.x)/deta) -
                                          (dgrid_ij.y_der.x/dxi)) + (Dx_Dy_ij.x/dxi);
    detady                              = mult * ((dgrid_ij.x_der.x * dDydy) -
                                                  (dgrid_ij.x_der.y * dDxdy) -
                                                  (ONEHALF * temp));


    /* Initialize to 0 */
    ders                                = zero (bh_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = ONE;  // dF(x(i, j))/dx(i, j)
    ders.ddx.p                          = -dxidx + ((ONEHALF * c * n_ij.y)/deta);   // dF(p(i, j))/dx(i, j)
    ders.ddx.q                          = -detadx;  // dF(q(i, j))/dx(i, j)

    ders.ddy.y                          = ONE;  // dF(y(i, j))/dy(i, j)
    ders.ddy.p                          = -dxidy - ((ONEHALF * c * n_ij.x)/deta);   // dF(p(i, j))/dy(i, j)
    ders.ddy.q                          = -detady;  // dF(q(i, j))/dy(i, j)

    ders.ddp.p                          = ONE;  // dF(p(i, j))/dp(i, j)

    ders.ddq.q                          = ONE;  // dF(q(i, j))/dq(i, j)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side

   Input parameters: i1         - i index of point for which derivatives
                                  are being computed
                     j1         - j index of point for which derivatives
                                  are being computed
                     i2         - i index of point with respect to which
                                  derivatives are being computed
                     j2         - j index of point with respect to which
                                  derivatives are being computed
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     dgrid_ij   - values of first-order derivatives of
                                  (x, y) at point (i1, j1)
                     coeff_ij   - values of first-order derivatives
                                  coefficients at point (i1, j1)
                     d2grid_ij  - values of second-order derivatives of
                                  (x, y) at point (i1, j1)
                     Dx_Dy_ij   - Dx and Dy at point (i1, j1)
                     nablas_ij  - Laplacians of xi and eta at point (i1, j1)
                     n_ij       - outward normal vector at point (i1, j1)
*/
bh_jacobian_2D xi1_eta1_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                    long double deta, long double c, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij, point_2D Dx_Dy_ij,
                                    point_2D nablas_ij, point_2D n_ij)
{
    /* Return dFdu */
    bh_jacobian_2D      dFdu;


    if (i2 == i1 && j2 == j1 - 3)   // dF(u(i1, j1))/du(i1, j1 - 3)
    {
        dFdu                            = xi1_eta1_i_jm3 (deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 2 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 - 2, j1 - 2)
    {
        dFdu                            = xi1_eta1_im2_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1 - 1, j1 - 2)
    {
        dFdu                            = xi1_eta1_im1_jm2 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 2)    // dF(u(i1, j1))/du(i1, j1 - 2)
    {
        dFdu                            = xi1_eta1_i_jm2 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 2 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 - 2, j1 - 1)
    {
        dFdu                            = xi1_eta1_im2_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1 - 1, j1 - 1)
    {
        dFdu                            = xi1_eta1_im1_jm1 (dxi, deta, dgrid_ij, coeff_ij);
    } else if (i2 == i1 && j2 == j1 - 1)    // dF(u(i1, j1))/du(i1, j1 - 1)
    {
        dFdu                            = xi1_eta1_i_jm1 (dxi, deta, c, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else if (i2 == i1 - 3 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 3, j1)
    {
        dFdu                            = xi1_eta1_im3_j (dxi, dgrid_ij, coeff_ij);
    } else if (i2 == i1 - 2 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 2, j1)
    {
        dFdu                            = xi1_eta1_im2_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 - 1 && j2 == j1)    // dF(u(i1, j1))/du(i1 - 1, j1)
    {
        dFdu                            = xi1_eta1_im1_j (dxi, deta, dgrid_ij, coeff_ij,
                                                          d2grid_ij, Dx_Dy_ij, nablas_ij);
    } else if (i2 == i1 && j2 == j1)    // dF(u(i1, j1))/du(i1, j1)
    {
        dFdu                            = xi1_eta1_i_j (dxi, deta, c, dgrid_ij, coeff_ij,
                                                        d2grid_ij, Dx_Dy_ij, nablas_ij, n_ij);
    } else
    {
        dFdu                            = zero (bh_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for elliptic grid generation     */
/*     Equivalent to setting (p, q) as 0                     */
/*                                                           */
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------o------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*-----------------------------------------------------------*/
/*                Interior point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i - 1, j - 1) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     deta       - computational spacing in eta
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
ell_jacobian_2D ell_interior_im1_jm1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j - 1)
    ders.ddy.y                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i, j - 1) for elliptic grid
   generation

   Input parameters: deta       - computational spacing in eta
                     dgrid_ij   - first derivatives at point (i, j)
                     d2grid_ij  - second derivatives at point (i, j)
*/
ell_jacobian_2D ell_interior_i_jm1 (long double deta, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/deta;
    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.gamma/deta)); // dF(x(i, j))/dx(i, j - 1)
    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x)); // dF(x(i, j))/dy(i, j - 1)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y)); // dF(y(i, j))/dx(i, j - 1)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.gamma/deta)); // dF(y(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i + 1, j - 1) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     deta       - computational spacing in eta
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
ell_jacobian_2D ell_interior_ip1_jm1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i + 1, j - 1)
    ders.ddy.y                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i + 1, j - 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i - 1, j) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     dgrid_ij   - first derivatives at point (i, j)
                     d2grid_ij  - second derivatives at point (i, j)
*/
ell_jacobian_2D ell_interior_im1_j (long double dxi, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/dxi;
    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (coeff_ij.alpha/dxi)); // dF(x(i, j))/dx(i - 1, j)
    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.x)); // dF(x(i, j))/dy(i - 1, j)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.y)); // dF(y(i, j))/dx(i - 1, j)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (coeff_ij.alpha/dxi)); // dF(y(i, j))/dy(i - 1, j)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i, j) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     deta       - computational spacing in eta
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
ell_jacobian_2D ell_interior_i_j (long double dxi, long double deta,
                                  coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(x(i, j))/dx(i, j)
    ders.ddy.y                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(y(i, j))/dy(i, j)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i + 1, j) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     dgrid_ij   - first derivatives at point (i, j)
                     d2grid_ij  - second derivatives at point (i, j)
*/
ell_jacobian_2D ell_interior_ip1_j (long double dxi, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/dxi;
    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (coeff_ij.alpha/dxi)); // dF(x(i, j))/dx(i + 1, j)
    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.x)); // dF(x(i, j))/dy(i + 1, j)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.y)); // dF(y(i, j))/dx(i + 1, j)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (coeff_ij.alpha/dxi)); // dF(y(i, j))/dy(i + 1, j)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i - 1, j + 1) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     deta       - computational spacing in eta
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
ell_jacobian_2D ell_interior_im1_jp1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j + 1)
    ders.ddy.y                          = (coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j + 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i, j + 1) for elliptic grid
   generation

   Input parameters: deta       - computational spacing in eta
                     dgrid_ij   - first derivatives at point (i, j)
                     d2grid_ij  - second derivatives at point (i, j)
*/
ell_jacobian_2D ell_interior_i_jp1 (long double deta, grid_der_2D dgrid_ij,
                                    coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/deta;
    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.gamma/deta)); // dF(x(i, j))/dx(i, j + 1)
    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x)); // dF(x(i, j))/dy(i, j + 1)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y)); // dF(y(i, j))/dx(i, j + 1)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.gamma/deta)); // dF(y(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point (i + 1, j + 1) for elliptic grid
   generation

   Input parameters: dxi        - computational spacing in xi
                     deta       - computational spacing in eta
                     coeff_ij   - first derivative coefficients at point (i, j)
*/
ell_jacobian_2D ell_interior_ip1_jp1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(x(i, j))/dx(i + 1, j + 1)
    ders.ddy.y                          = -(coeff_ij.beta)/(TWO * dxi * deta); // dF(y(i, j))/dy(i + 1, j + 1)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for xi = 1
   points which constitute the matrix
   for Newton's method on the left hand
   side
*/
ell_jacobian_2D ell_interior_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                         long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                         grid_dder_2D d2grid_ij)
{
    /* Return dFdu */
    ell_jacobian_2D         dFdu;


    if (i2 == i1 - 1 && j2 == j1 - 1)
    {
        dFdu                            = ell_interior_im1_jm1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 && j2 == j1 - 1)
    {
        dFdu                            = ell_interior_i_jm1 (deta, dgrid_ij, coeff_ij, d2grid_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1 - 1)
    {
        dFdu                            = ell_interior_ip1_jm1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 - 1 && j2 == j1)
    {
        dFdu                            = ell_interior_im1_j (dxi, dgrid_ij, coeff_ij, d2grid_ij);
    }
    else if (i2 == i1 && j2 == j1)
    {
        dFdu                            = ell_interior_i_j (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1)
    {
        dFdu                            = ell_interior_ip1_j (dxi, dgrid_ij, coeff_ij, d2grid_ij);
    }
    else if (i2 == i1 - 1 && j2 == j1 + 1)
    {
        dFdu                            = ell_interior_im1_jp1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 && j2 == j1 + 1)
    {
        dFdu                            = ell_interior_i_jp1 (deta, dgrid_ij, coeff_ij, d2grid_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1 + 1)
    {
        dFdu                            = ell_interior_ip1_jp1 (dxi, deta, coeff_ij);
    }
    else
    {
        dFdu                            = zero (ell_jacobian_2D);
    }


    return dFdu;
}






/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for elliptic grid generation     */
/*     with non-zero control functions (p, q)                */
/*                                                           */
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/*                      x------x------x                      */
/*                      |      |      |                      */
/*                      x------o------x                      */
/*                      |      |      |                      */
/*                      x------x------x                      */
/*-----------------------------------------------------------*/
/*                Interior point derivatives                 */
/*-----------------------------------------------------------*/



/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point u(i - 1, j - 1) for Poisson
   grid generation

   Input parameters: dxi      - spacing of xi in computational grid
                     deta     - spacing of eta in computational grid
                     coeff_ij - first derivative coefficients at
                                point (i, j)
*/
ell_jacobian_2D psn_interior_im1_jm1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -coeff_ij.beta/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j - 1)
    ders.ddy.y                          = -coeff_ij.beta/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
  Compute derivatives of F(u(i, j)) with respect
  to interior point u(i, j - 1) for Poisson grid
  generation

  Input parameters: deta        - spacing of eta in computational grid
                    dgrid_ij    - derivatives of (x, y) with respect to (xi, eta)
                                  at point (i, j)
                    coeff_ij    - first derivative coefficients at point (i, j)
                    d2grid_ij   - second-order derivatives of (x, y) with respect to
                                  (xi, eta) at point (i, j)
                    pq_ij       - values of control functions (p, q) at point (i, j)
*/
ell_jacobian_2D psn_interior_i_jm1 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                    grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult, tempx, tempy, temp;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/deta;
    tempx                               = (pq_ij.x * dgrid_ij.x_der.x) + (pq_ij.y * dgrid_ij.y_der.x);
    tempy                               = (pq_ij.x * dgrid_ij.x_der.y) + (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.gamma/deta) - ((coeff_ij.J * coeff_ij.J * pq_ij.y)/TWO);

    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * tempx) + temp); // dF(x(i, j))/dx(i, j - 1)
    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * tempx)); // dF(x(i, j))/dy(i, j - 1)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.y * tempy)); // dF(y(i, j))/dx(i, j - 1)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.x_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.x * tempy) + temp); // dF(y(i, j))/dy(i, j - 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point u(i + 1, j - 1) for Poisson
   grid generation

   Input parameters: dxi      - spacing of xi in computational grid
                     deta     - spacing of eta in computational grid
                     coeff_ij - first derivative coefficients at
                                point (i, j)
*/
ell_jacobian_2D psn_interior_ip1_jm1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = coeff_ij.beta/(TWO * dxi * deta); // dF(x(i, j))/dx(i + 1, j - 1)
    ders.ddy.y                          = coeff_ij.beta/(TWO * dxi * deta); // dF(y(i, j))/dy(i + 1, j - 1)


    return ders;
}






/*
  Compute derivatives of F(u(i, j)) with respect
  to interior point u(i - 1, j) for Poisson grid
  generation

  Input parameters: deta        - spacing of eta in computational grid
                    dgrid_ij    - derivatives of (x, y) with respect to (xi, eta)
                                  at point (i, j)
                    coeff_ij    - first derivative coefficients at point (i, j)
                    d2grid_ij   - second-order derivatives of (x, y) with respect to
                                  (xi, eta) at point (i, j)
                    pq_ij       - values of control functions (p, q) at point (i, j)
*/
ell_jacobian_2D psn_interior_im1_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                    grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult, tempx, tempy, temp;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/dxi;
    tempx                               = (pq_ij.x * dgrid_ij.x_der.x) + (pq_ij.y * dgrid_ij.y_der.x);
    tempy                               = (pq_ij.x * dgrid_ij.x_der.y) + (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.alpha/dxi) - ((coeff_ij.J * coeff_ij.J * pq_ij.x)/TWO);

    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * tempx) + temp); // dF(x(i, j))/dx(i - 1, j)
    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * tempx)); // dF(x(i, j))/dy(i - 1, j)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.xy_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.y * tempy)); // dF(y(i, j))/dx(i - 1, j)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.xy_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.x * tempy) + temp); // dF(y(i, j))/dy(i - 1, j)

    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point u(i, j) for Poisson
   grid generation

   Input parameters: dxi      - spacing of xi in computational grid
                     deta     - spacing of eta in computational grid
                     coeff_ij - first derivative coefficients at
                                point (i, j)
*/
ell_jacobian_2D psn_interior_i_j (long double dxi, long double deta, coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(x(i, j))/dx(i, j)
    ders.ddy.y                          = -TWO * ((coeff_ij.alpha/(dxi * dxi)) +
                                                  (coeff_ij.gamma/(deta * deta))); // dF(y(i, j))/dy(i, j)


    return ders;
}






/*
  Compute derivatives of F(u(i, j)) with respect
  to interior point u(i + 1, j) for Poisson grid
  generation

  Input parameters: deta        - spacing of eta in computational grid
                    dgrid_ij    - derivatives of (x, y) with respect to (xi, eta)
                                  at point (i, j)
                    coeff_ij    - first derivative coefficients at point (i, j)
                    d2grid_ij   - second-order derivatives of (x, y) with respect to
                                  (xi, eta) at point (i, j)
                    pq_ij       - values of control functions (p, q) at point (i, j)
*/
ell_jacobian_2D psn_interior_ip1_j (long double dxi, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                    grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult, tempx, tempy, temp;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/dxi;
    tempx                               = (pq_ij.x * dgrid_ij.x_der.x) + (pq_ij.y * dgrid_ij.y_der.x);
    tempy                               = (pq_ij.x * dgrid_ij.x_der.y) + (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.alpha/dxi) + ((coeff_ij.J * coeff_ij.J * pq_ij.x)/TWO);

    ders.ddx.x                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.x) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.x) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * tempx) + temp); // dF(x(i, j))/dx(i + 1, j)
    ders.ddy.x                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.x) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.x) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * tempx)); // dF(x(i, j))/dy(i + 1, j)
    ders.ddx.y                          = mult * (-(dgrid_ij.y_der.x * d2grid_ij.xy_der.y) +
                                                  (dgrid_ij.x_der.x * d2grid_ij.y_der.y) +
                                                  (coeff_ij.J * dgrid_ij.y_der.y * tempy)); // dF(y(i, j))/dx(i + 1, j)
    ders.ddy.y                          = mult * (-(dgrid_ij.y_der.y * d2grid_ij.xy_der.y) +
                                                  (dgrid_ij.x_der.y * d2grid_ij.y_der.y) -
                                                  (coeff_ij.J * dgrid_ij.y_der.x * tempy) + temp); // dF(y(i, j))/dy(i + 1, j)

    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point u(i - 1, j + 1) for Poisson
   grid generation

   Input parameters: dxi      - spacing of xi in computational grid
                     deta     - spacing of eta in computational grid
                     coeff_ij - first derivative coefficients at
                                point (i, j)
*/
ell_jacobian_2D psn_interior_im1_jp1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = coeff_ij.beta/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j + 1)
    ders.ddy.y                          = coeff_ij.beta/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j + 1)


    return ders;
}






/*
  Compute derivatives of F(u(i, j)) with respect
  to interior point u(i, j + 1) for Poisson grid
  generation

  Input parameters: deta        - spacing of eta in computational grid
                    dgrid_ij    - derivatives of (x, y) with respect to (xi, eta)
                                  at point (i, j)
                    coeff_ij    - first derivative coefficients at point (i, j)
                    d2grid_ij   - second-order derivatives of (x, y) with respect to
                                  (xi, eta) at point (i, j)
                    pq_ij       - values of control functions (p, q) at point (i, j)
*/
ell_jacobian_2D psn_interior_i_jp1 (long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                    grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;

    /* Local variables */
    long double             mult, tempx, tempy, temp;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    mult                                = ONE/deta;
    tempx                               = (pq_ij.x * dgrid_ij.x_der.x) + (pq_ij.y * dgrid_ij.y_der.x);
    tempy                               = (pq_ij.x * dgrid_ij.x_der.y) + (pq_ij.y * dgrid_ij.y_der.y);
    temp                                = (coeff_ij.gamma/deta) + ((coeff_ij.J * coeff_ij.J * pq_ij.y)/TWO);

    ders.ddx.x                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.x) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * tempx) + temp); // dF(x(i, j))/dx(i, j + 1)
    ders.ddy.x                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.x) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.x) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * tempx)); // dF(x(i, j))/dy(i, j + 1)
    ders.ddx.y                          = mult * ((dgrid_ij.y_der.x * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.x * d2grid_ij.xy_der.y) -
                                                  (coeff_ij.J * dgrid_ij.x_der.y * tempy)); // dF(y(i, j))/dx(i, j + 1)
    ders.ddy.y                          = mult * ((dgrid_ij.y_der.y * d2grid_ij.x_der.y) -
                                                  (dgrid_ij.x_der.y * d2grid_ij.xy_der.y) +
                                                  (coeff_ij.J * dgrid_ij.x_der.x * tempy) + temp); // dF(y(i, j))/dy(i, j + 1)


    return ders;
}






/*
   Compute derivatives of F(u(i, j)) with respect
   to interior point u(i + 1, j + 1) for Poisson
   grid generation

   Input parameters: dxi      - spacing of xi in computational grid
                     deta     - spacing of eta in computational grid
                     coeff_ij - first derivative coefficients at
                                point (i, j)
*/
ell_jacobian_2D psn_interior_ip1_jp1 (long double dxi, long double deta,
                                      coeffs_1 coeff_ij)
{
    /* Return ders */
    ell_jacobian_2D         ders;


    /* Initialize to 0 */
    ders                                = zero (ell_jacobian_2D);

    /* Non-zero derivatives */
    ders.ddx.x                          = -coeff_ij.beta/(TWO * dxi * deta); // dF(x(i, j))/dx(i - 1, j - 1)
    ders.ddy.y                          = -coeff_ij.beta/(TWO * dxi * deta); // dF(y(i, j))/dy(i - 1, j - 1)


    return ders;
}






/*
   Function to compute the derivatives
   dF(u(i1, j1))/du(i2, j2) for interior
   points which constitute the matrix
   for Newton's method on the left hand
   side
*/
ell_jacobian_2D psn_interior_point_ders (int i1, int j1, int i2, int j2, long double dxi,
                                         long double deta, grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                         grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return dFdu */
    ell_jacobian_2D         dFdu;


    if (i2 == i1 - 1 && j2 == j1 - 1)
    {
        dFdu                            = psn_interior_im1_jm1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 && j2 == j1 - 1)
    {
        dFdu                            = psn_interior_i_jm1 (deta, dgrid_ij, coeff_ij, d2grid_ij, pq_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1 - 1)
    {
        dFdu                            = psn_interior_ip1_jm1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 - 1 && j2 == j1)
    {
        dFdu                            = psn_interior_im1_j (dxi, dgrid_ij, coeff_ij, d2grid_ij, pq_ij);
    }
    else if (i2 == i1 && j2 == j1)
    {
        dFdu                            = psn_interior_i_j (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1)
    {
        dFdu                            = psn_interior_ip1_j (dxi, dgrid_ij, coeff_ij, d2grid_ij, pq_ij);
    }
    else if (i2 == i1 - 1 && j2 == j1 + 1)
    {
        dFdu                            = psn_interior_im1_jp1 (dxi, deta, coeff_ij);
    }
    else if (i2 == i1 && j2 == j1 + 1)
    {
        dFdu                            = psn_interior_i_jp1 (deta, dgrid_ij, coeff_ij, d2grid_ij, pq_ij);
    }
    else if (i2 == i1 + 1 && j2 == j1 + 1)
    {
        dFdu                            = psn_interior_ip1_jp1 (dxi, deta, coeff_ij);
    }
    else
    {
        dFdu                            = zero (ell_jacobian_2D);
    }


    return dFdu;
}
