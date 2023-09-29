*NAME:         bm23
*ROWS:         20
*COLUMNS:      27
*INTEGER:      27
*NONZERO:      478
*BEST SOLN:    34 (opt)
*LP SOLN:      20.57
*SOURCE:       B. Bouvier and G. Messoumain (Univ. of Grenoble) 
*              Harlan Crowder (IBM)
*              E. Andrew Boyd (Rice University)
*APPLICATION:  unknown
*COMMENTS:     pure 0/1 IP                      
*      
*      
NAME          BM23
ROWS
 N  R100    
 L  R101    
 L  R102    
 L  R103    
 L  R104    
 L  R105    
 L  R106    
 L  R107    
 L  R108    
 L  R109    
 L  R110    
 L  R111    
 L  R112    
 L  R113    
 L  R114    
 L  R115    
 L  R116    
 L  R117    
 L  R118    
 L  R119    
 L  R120    
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    C101      R100                 3   R101                 2
    C101      R102                -6   R103                 5
    C101      R104                -9   R105                 8
    C101      R106                 7   R107                -1
    C101      R108                -1   R109                 4
    C101      R111                -1   R112                 8
    C101      R113                -2   R114                 7
    C101      R115                 7   R116                 3
    C101      R117                -3   R118                 5
    C101      R119                -4   R120                -6
    C102      R100                 2   R101                 6
    C102      R102                 5   R103                -6
    C102      R104                -5   R105                -7
    C102      R106                 5   R107                 3
    C102      R108                 2   R109                -6
    C102      R111                 3   R113                -7
    C102      R114                 9   R115                 1
    C102      R116                 6   R117                -1
    C102      R118                -9   R119                 4
    C102      R120                 8
    C103      R100                 5   R101                 1
    C103      R102                -8   R103                -5
    C103      R105                 4   R106                 2
    C103      R107                 3   R108                -6
    C103      R109                -7   R111                 3
    C103      R112                 5   R113                -8
    C103      R114                -7   R115                 4
    C103      R116                 1   R117                 6
    C103      R118                -5   R119                -8
    C103      R120                 2
    C104      R100                 8   R102                -3
    C104      R103                -3   R104                 9
    C104      R105                 5   R107                 4
    C104      R108                -9   R109                -2
    C104      R111                -4   R113                 6
    C104      R114                -5   R115                -3
    C104      R116                -3   R117                -6
    C104      R119                -5
    C105      R100                 6   R101                 3
    C105      R103                -8   R104                -1
    C105      R105                 9   R106                 6
    C105      R107                 1   R109                -2
    C105      R111                 5   R112                 2
    C105      R113                -2   R114                 1
    C105      R116                -5   R117                 1
    C105      R119                -7   R120                -8
    C106      R100                 9   R101                 3
    C106      R102                -1   R103                 8
    C106      R104                 8   R105                -1
    C106      R106                 6   R108                 7
    C106      R111                 5   R112                 6
    C106      R113                -5   R114                -5
    C106      R115                -6   R116                 6
    C106      R117                -2   R118                -6
    C106      R119                -7   R120                -7
    C107      R100                 7   R101                 2
    C107      R102                -3   R103                -9
    C107      R104                -3   R105                 7
    C107      R106                 7   R107                 4
    C107      R108                -9   R109                -6
    C107      R111                 7   R112                 7
    C107      R113                -9   R114                -5
    C107      R115                -3   R116                -9
    C107      R118                 7   R119                -4
    C107      R120                -6
    C108      R100                 4   R101                 6
    C108      R102                -8   R103                -2
    C108      R104                 9   R105                 1
    C108      R106                 6   R107                -1
    C108      R108                 9   R109                -6
    C108      R111                 8   R112                 1
    C108      R113                -2   R114                 4
    C108      R115                 7   R116                 6
    C108      R117                 3   R118                 1
    C108      R119                 5
    C109      R100                 5   R101                 2
    C109      R102                -9   R104                 9
    C109      R105                -3   R106                -7
    C109      R107                -6   R108                 6
    C109      R109                 7   R111                -8
    C109      R112                -1   R113                 2
    C109      R114                 9   R115                 8
    C109      R116                 3   R117                -3
    C109      R118                -1   R120                 9
    C110      R100                 6   R101                 2
    C110      R102                 3   R103                 9
    C110      R104                 3   R105                 2
    C110      R106                -7   R108                -4
    C110      R109                -4   R110                -1
    C110      R112                 1   R113                -8
    C110      R114                 3   R115                 1
    C110      R116                -9   R117                 3
    C110      R118                 7   R119                -8
    C110      R120                 4
    C111      R100                 3   R101                 5
    C111      R102                -8   R103                -1
    C111      R106                -1   R107                 8
    C111      R108                 5   R110                -1
    C111      R112                 4   R113                 1
    C111      R114                -4   R115                -6
    C111      R116                -6   R117                 6
    C111      R118                 2   R119                 9
    C111      R120                -8
    C112      R100                 2   R101                 3
    C112      R102                 6   R103                 7
    C112      R104                -7   R105                -3
    C112      R106                -8   R108                 3
    C112      R109                 2   R110                -1
    C112      R111                -1   R112                -5
    C112      R113                -3   R115                 6
    C112      R116                -3   R117                -1
    C112      R118                -6   R119                 9
    C113      R100                 4   R101                 3
    C113      R102                -3   R103                 9
    C113      R104                 5   R105                -5
    C113      R106                -3   R107                 1
    C113      R108                 1   R109                 8
    C113      R110                -1   R111                -2
    C113      R112                -7   R113                 5
    C113      R115                 8   R116                -6
    C113      R118                 2   R120                -3
    C114      R100                 5   R101                 7
    C114      R102                -8   R103                -9
    C114      R104                 4   R105                 9
    C114      R106                 9   R107                -5
    C114      R108                -3   R110                 1
    C114      R111                 4   R112                -8
    C114      R113                 6   R114                 7
    C114      R115                -3   R116                -6
    C114      R117                 6   R118                -2
    C114      R119                 8   R120                 3
    C115      R100                 8   R101                 9
    C115      R102                -6   R103                 4
    C115      R104                -9   R105                 7
    C115      R106                 1   R107                -4
    C115      R108                -9   R109                 4
    C115      R110                 1   R111                -6
    C115      R112                 2   R113                -8
    C115      R114                -7   R115                 8
    C115      R116                 6   R117                 6
    C115      R118                -4   R119                -8
    C116      R100                 7   R101                 4
    C116      R102                 7   R103                 2
    C116      R104                -1   R105                 2
    C116      R106                -6   R107                -1
    C116      R108                -1   R109                -4
    C116      R110                -1   R111                 6
    C116      R112                -8   R113                 1
    C116      R114                -9   R115                -5
    C116      R116                 2   R117                 8
    C116      R118                -5   R119                -4
    C116      R120                 4
    C117      R100                 3   R101                 2
    C117      R102                 6   R103                -2
    C117      R104                -9   R105                 5
    C117      R106                 6   R107                -9
    C117      R108                -6   R109                 3
    C117      R110                 1   R111                 7
    C117      R112                -5   R113                 1
    C117      R115                 1   R116                -7
    C117      R117                -2   R118                 4
    C117      R119                -8   R120                -6
    C118      R100                 5   R101                 2
    C118      R102                -2   R103                -9
    C118      R104                 7   R105                 4
    C118      R106                 6   R107                -7
    C118      R108                -8   R109                 2
    C118      R110                 1   R111                -8
    C118      R112                -5   R113                -9
    C118      R114                 9   R115                 8
    C118      R116                -6   R117                -4
    C118      R118                 5   R119                 7
    C118      R120                -2
    C119      R100                 2   R101                 3
    C119      R102                 3   R104                -3
    C119      R105                -7   R106                 1
    C119      R107                 2   R108                -9
    C119      R109                -3   R110                 1
    C119      R111                 8   R112                 9
    C119      R113                 8   R114                -3
    C119      R115                 4   R117                 1
    C119      R118                -1   R119                 4
    C119      R120                 3
    C120      R100                 4   R101                 4
    C120      R102                -7   R103                -3
    C120      R104                -7   R105                -9
    C120      R106                 1   R107                 2
    C120      R108                -9   R110                 2
    C120      R113                 3   R115                -8
    C120      R116                -7   R118                -8
    C120      R119                -6   R120                 5
    C121      R100                 3   R101                 9
    C121      R102                 4   R103                 3
    C121      R104                 4   R105                 3
    C121      R106                 6   R107                 1
    C121      R108                 8   R109                 4
    C121      R110                 2   R111                -1
    C121      R112                -1   R114                -9
    C121      R115                -3   R116                -4
    C121      R118                -2   R119                 9
    C121      R120                 3
    C122      R100                 2   R101                 9
    C122      R102                 5   R104                 1
    C122      R105                 2   R106                -3
    C122      R107                -3   R108                 5
    C122      R109                 3   R110                 2
    C122      R111                -2   R112                -9
    C122      R113                -8   R115                 8
    C122      R116                -6   R117                 6
    C122      R118                -2   R119                 5
    C122      R120                -4
    C123      R100                 5   R101                 9
    C123      R102                 5   R103                -4
    C123      R104                -1   R105                -6
    C123      R106                 5   R107                 7
    C123      R109                 7   R110                -2
    C123      R111                -4   R112                -3
    C123      R113                -3   R114                 9
    C123      R115                -1   R116                 5
    C123      R117                 1   R118                -4
    C123      R119                -7   R120                -9
    C124      R100                 6   R101                 2
    C124      R102                 5   R103                 7
    C124      R104                 1   R105                 8
    C124      R107                 8   R108                 5
    C124      R109                 9   R110                 2
    C124      R111                -7   R112                -6
    C124      R113                -5   R114                -6
    C124      R116                -9   R117                 1
    C124      R118                 1   R119                 3
    C124      R120                 9
    C125      R100                 8   R101                 8
    C125      R102                -3   R103                -1
    C125      R104                 3   R105                 6
    C125      R106                 5   R107                -5
    C125      R108                 8   R109                 8
    C125      R110                -2   R111                -7
    C125      R112                -1   R113                -2
    C125      R115                 5   R116                 5
    C125      R117                 9   R119                -6
    C125      R120                -4
    C126      R100                 1   R101                 6
    C126      R102                 3   R103                 4
    C126      R105                -8   R106                 6
    C126      R107                -2   R108                -3
    C126      R109                -3   R110                -2
    C126      R111                 7   R112                 5
    C126      R113                -2   R114                -4
    C126      R116                 4   R117                -7
    C126      R118                -2   R120                -5
    C127      R100                 1   R101                 5
    C127      R103                -6   R105                -9
    C127      R106                -8   R107                -1
    C127      R108                -1   R109                 8
    C127      R110                 2   R111                 8
    C127      R112                -4   R113                 8
    C127      R114                 9   R115                 2
    C127      R116                 2   R117                 1
    C127      R118                 1   R119                -5
    C127      R120                 1
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS       R101                65   R102               -15
    RHS       R103               -10   R104                32
    RHS       R105                14   R106                33
    RHS       R107                -5   R108                 3
    RHS       R109                18   R110                10
    RHS       R111                17   R112                -5
    RHS       R113                 4   R114                21
    RHS       R115                24   R116                -6
    RHS       R117                17   R118               -13
    RHS       R119               -30   R120                 1
BOUNDS
 UP ONE       C101                 1
 UP ONE       C102                 1
 UP ONE       C103                 1
 UP ONE       C104                 1
 UP ONE       C105                 1
 UP ONE       C106                 1
 UP ONE       C107                 1
 UP ONE       C108                 1
 UP ONE       C109                 1
 UP ONE       C110                 1
 UP ONE       C111                 1
 UP ONE       C112                 1
 UP ONE       C113                 1
 UP ONE       C114                 1
 UP ONE       C115                 1
 UP ONE       C116                 1
 UP ONE       C117                 1
 UP ONE       C118                 1
 UP ONE       C119                 1
 UP ONE       C120                 1
 UP ONE       C121                 1
 UP ONE       C122                 1
 UP ONE       C123                 1
 UP ONE       C124                 1
 UP ONE       C125                 1
 UP ONE       C126                 1
 UP ONE       C127                 1
ENDATA
