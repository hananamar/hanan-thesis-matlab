function res = catalogueK(a,b,x,y)
    D = det([ 1 x(1) y(1) ;...
              1 x(2) y(2) ;...
              1 x(3) y(3) ]);
          
    A = [ 1 y(2)-y(3) x(3)-x(2) ;...
          1 y(3)-y(1) x(1)-x(3) ;...
          1 y(1)-y(2) x(2)-x(1) ]/D;
    
    res = D*(A(a,2)*A(b,2) + A(a,3)*A(b,3))/2;
end