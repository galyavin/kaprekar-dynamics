# Fixed box test and predictions

For initial difference state q=(a,b,c,d,e), use the normalized chart
x=(a,b,c,2(b-e)-B,B-c-d)/B and clock t=T_B/B,
let theta(x)=121049/1000000 + sum_i v_i (x_i-u_i)/s_i. The integer arrays below
are in millionths, so each entry is divided by 1000000 before use.

- u = (940800, 738087, 579763, 484, 2491)
- s = (446, 365, 342, 280, 449)
- v = (590, -307, -447, 127, -296)

All coefficients are exact. The rational coefficients were fixed before route verification.
The fixed test is F(x,t)=min(1,max(0,500*(t-theta(x)+1/1000))). The chart is
linear in the initial gaps, so this defines a bounded Lipschitz test of the
joint observation, including away from the box.

Predicted whole-box claims for each synchronized bucket, once B is large:

- The residual route clock groups in increasing order are S, M1, M2 and L.
- For B=12j+5 their counts are (542,52,82,380) and for B+6 they are
  (598,52,94,312), out of 1056 per arm. The contrast order is B minus B+6.
- Each S, M1 and M2 class gives F=0. Each L class gives F=1.
- The limiting contrast under the two arms' own conditional word weights
  is 68/28512 = 17/7128, conditional on the whole-box guard certificate.

The test and predictions were fixed before this certificate and its replay.
