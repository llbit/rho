expected <- eval(parse(text="structure(\"with 5\\\\% of the range added to each end.\\n\", Rd_tag = \"TEXT\")"));                   
test(id=0, code={                   
argv <- eval(parse(text="list(\"([&$%_#])\", \"\\\\\\\\\\\\1\", structure(\"with 5% of the range added to each end.\\n\", Rd_tag = \"TEXT\"), FALSE, TRUE, FALSE, TRUE)"));                   
.Internal(gsub(argv[[1]], argv[[2]], argv[[3]], argv[[4]], argv[[5]], argv[[6]], argv[[7]]));                   
}, o=expected);                   

