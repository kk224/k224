# 224ori.txt
The system generated from Equation (5) and (6) is contained in '224ori.txt'. There are 576 polynomials (should equal to 0) in 512 + 64 = 576 variables. x0\~x511 stand for the unknowns in (a) of the 1st round in Figure 7. For example, A_{x=0, y=0, z=3}, A_{x=2, y=0, z=0}, and A_{x=0, y=1, z=0} are represented by x3, x64, and x128 respectively. x600\~x663 refer to the selected pair (p_{0, 3}, p_{0, 4}) in Equation (5) and (6).

Please note that Equation (10) is not included in the above system. We verify whether it holds after a first message block is achieved.


# 224param.txt
There are 448 linear polynomials and 128 quadratic polynomials in '224ori.txt'. Reducing these 128 quadratic polynomials by all linear polynomials leads to 128 new nonlinear polynomials, which are included in '224param.txt'. The reduction can be done by running './F4_param_b1' directly, or using any computer algebra system, e.g. Magma. The monomial ordering is x0 > x1 > .... So after the reduction, the variables x0~x447 are eliminated.



# 224rename.txt
For convenience, the variables x600-x663 in polynomials of '224param.txt' are renamed as x0~x63.
