expected <- eval(parse(text="c(69.0775527898214, 138.155105579643, 207.232658369464, 276.310211159285, 345.387763949107, 414.465316738928, 483.54286952875, 552.620422318571, 621.697975108392, 690.775527898214)"));         
test(id=0, code={         
argv <- eval(parse(text="list(c(1e+30, 1e+60, 1e+90, 1e+120, 1e+150, 1e+180, 1e+210, 1e+240, 1e+270, 1e+300))"));         
do.call(`digamma`, argv);         
}, o=expected);         

