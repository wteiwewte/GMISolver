NAME          TESTPROB
ROWS
 N  COST
 L  LIM1
 L  LIM2
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    XONE      LIM1                 3   LIM2                -3
    XTWO      LIM1                 2   LIM2                 2
    XTWO      COST                -1
    MARK0001  'MARKER'                 'INTEND'
RHS
    RHS1      LIM1                 6   LIM2                 0
BOUNDS
 UP BND1      XONE                 1000
 UP BND1      XTWO                 1000
ENDATA
