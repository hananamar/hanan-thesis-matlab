function res = catalogueK()
    global h;

    A = [ 1 -1;...
         -1  1];
     
    res = A/h;
end