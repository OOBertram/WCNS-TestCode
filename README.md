# WCNS-TestCode
## Introduction

This code includes a test suite for the 5th-order Weighted Compact Nonlinear Scheme (WCNS), which allows testing the accuracy and resolution of the 5th-order WCNS. Please cite the following articles when using this code:

1.Xuan Liu, Meiyuan Zhen, Jinsheng Cai, Fei Liao; A fast hybrid weighted compact nonlinear scheme with enhanced smoothness indicator. Physics of Fluids 1 March 2025; 37 (3): 036148. https://doi.org/10.1063/5.0259781

2.Liu, X., Min, Y., Cai, J., Ma, Y. and Yan, Z.-G. (2025), A Filtered Embedded Weighted Compact Non-Linear Scheme for Hyperbolic Conservation Law. Int J Numer Meth Fluids. https://doi.org/10.1002/fld.5366

 **Efficiency-Double Mach Refelection** 
| Case | Scheme    | CFL | Grid         | Time Use(s) | Time Ratio(%) | Shock Detect Time(s) | Detect Time Ratio (DTR)(%) | DTR-TOT |
|------|-----------|-----|--------------|-------------|---------------|-----------------------|-----------------------------|---------|
| DMR  | WCNS-Z    | 0.1 | 960×240      | 22393       | -             | -                     | -                           | -       |
| DMR  | WCNS-T    | 0.1 | 960×240      | 28193       | 125.9         | -                     | -                           | -       |
| DMR  | WCNS-CU6  | 0.1 | 960×240      | 26975       | 120.5         | -                     | -                           | -       |
| DMR  | WCNS-TOT  | 0.1 | 960×240      | 20044       | 89.51         | 16019                | 79.92                       | 1.000   |
| DMR  | WCNS-Ada   | 0.1 | 960×240      | 31328       | 139.9         | 6218                 | 19.85                       | 2.576   |
| DMR  | WCNS-FT   | 0.1 | 960×240      | 4059        | 18.13         | 32.7                 | 0.8050                      | 489.9   |
| DMR  | WCNS-Z    | 0.2 | 960×240      | 11086       | -             | -                     | -                           | -       |
| DMR  | WCNS-T    | 0.2 | 960×240      | 13957       | 125.9         | -                     | -                           | -       |
| DMR  | WCNS-CU6  | 0.2 | 960×240      | 13358       | 120.5         | -                     | -                           | -       |
| DMR  | WCNS-TOT  | 0.2 | 960×240      | 9972        | 89.95         | 8009.9               | 80.32                       | 1.000   |
| DMR  | WCNS-Ada   | 0.2 | 960×240      | 15542       | 140.2         | 3108                 | 20.00                       | 2.577   |
| DMR  | WCNS-FT   | 0.2 | 960×240      | 2066        | 18.64         | 32.0                 | 1.550                       | 250.3   |
| DMR  | WCNS-Z    | 0.3 | 960×240      | 7294        | -             | -                     | -                           | -       |
| DMR  | WCNS-T    | 0.3 | 960×240      | 9176        | 125.8         | -                     | -                           | -       |
| DMR  | WCNS-CU6  | 0.3 | 960×240      | 8755        | 120.0         | -                     | -                           | -       |
| DMR  | WCNS-TOT  | 0.3 | 960×240      | 6628        | 90.87         | 5317.5               | 80.22                       | 1.000   |
| DMR  | WCNS-Ada   | 0.3 | 960×240      | 10204       | 139.9         | 2092                 | 20.502                      | 2.542   |
| DMR  | WCNS-FT   | 0.3 | 960×240      | 1402        | 19.22         | 31.9                 | 2.272                       | 166.7   |
| DMR  | WCNS-Z    | 0.4 | 960×240      | 5457        | -             | -                     | -                           | -       |
| DMR  | WCNS-T    | 0.4 | 960×240      | 6876        | 126.0         | -                     | -                           | -       |
| DMR  | WCNS-CU6  | 0.4 | 960×240      | 6595        | 120.9         | -                     | -                           | -       |
| DMR  | WCNS-TOT  | 0.4 | 960×240      | 4945        | 90.62         | 3965.7               | 80.18                       | 1.000   |
| DMR  | WCNS-Ada   | 0.4 | 960×240      | 7684        | 140.8         | 1477                 | 19.22                       | 2.685   |
| DMR  | WCNS-FT   | 0.4 | 960×240      | 1072        | 19.64         | 31.5                 | 2.940                       | 125.9   |
| DMR  | WCNS-Z    | 0.5 | 960×240      | 4376        | -             | -                     | -                           | -       |
| DMR  | WCNS-T    | 0.5 | 960×240      | 5438        | 125.2         | -                     | -                           | -       |
| DMR  | WCNS-CU6  | 0.5 | 960×240      | 5261        | 120.2         | -                     | -                           | -       |
| DMR  | WCNS-TOT  | 0.5 | 960×240      | 3934        | 89.56         | 3145                 | 79.94                       | 1.000   |
| DMR  | WCNS-Ada   | 0.5 | 960×240      | 6146        | 140.1         | 1182                 | 19.23                       | 2.661   |
| DMR  | WCNS-FT   | 0.5 | 960×240      | 885.5       | 19.77         | 31.3                 | 3.538                       | 100.5   |


 **CPU_INFO** 

    Intel(R) Core(TM) i7-14700KF
    CPU benchmark speed: 3.40 GHz

## Usage

gfortran *.f90 -o a.out -freal-4-real-8 -fno-align-commons

## Issues and Support

If you have any questions or encounter any bugs in the program, please contact us at: bertram0507@163.com or bertram0298@gmail.com.
