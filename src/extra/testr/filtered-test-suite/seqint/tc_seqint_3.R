expected <- eval(parse(text="1L"));            
test(id=0, code={            
argv <- eval(parse(text="list(FALSE)"));            
do.call(`seq.int`, argv);            
}, o=expected);            

