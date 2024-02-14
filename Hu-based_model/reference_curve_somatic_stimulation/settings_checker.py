def LS_AC__value_checker(vLeftShift, AC):
    if vLeftShift!=0.0 and AC==0.0:
        print('ERROR: AC=0 HENCE vLeftShift HAS NO EFFECT!!!')
        exit()
    if AC<0.0:
        print('ERROR: AC<0 is NONSENSICAL, FRACTION OF LEFT-SHIFTED CHANNELS (AC) MUST BE POSITIVE!!!')
        exit()

