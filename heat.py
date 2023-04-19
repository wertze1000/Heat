#(F -> function evaluated at t1 or t2 to find F at tEval)
def LinReg(Ft1, Ft2, t1, t2, tEval):
    slope = (Ft2 - Ft1) / (t2 - t1)
    p = Ft2 - slope*t2
    FtEval = slope*tEval + p
    return FtEval