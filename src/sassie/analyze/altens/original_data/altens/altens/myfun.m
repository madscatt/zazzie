function F=myfun(x,C);





ar=x(1,1);
br=x(1,2);
gr=x(1,3);


F=[(cos(ar)*cos(br)*cos(gr)-sin(ar)*sin(gr))-C(1,1);(sin(ar)*cos(br)*cos(gr)+cos(ar)*sin(gr))-C(1,2);...
        (-sin(br)*cos(gr))-C(1,3);(-cos(ar)*cos(br)*sin(gr)-sin(ar)*cos(gr))-C(2,1);(-sin(ar)*cos(br)*sin(gr)+cos(ar)*cos(gr))-C(2,2);...
        (sin(br)*sin(gr))-C(2,3);(cos(ar)*sin(br))-C(3,1);(sin(ar)*sin(br))-C(3,2);...
        cos(br)-C(3,3)];
    
    
return



