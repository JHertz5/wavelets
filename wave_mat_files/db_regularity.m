for i=1:15
    [phi,psi,xval]=wavefun('db3',i);
    stairs(xval,phi)
    xlabel(i)
    pause(1)
end
pause(5)
figure
for i=1:15
    [phi,psi,xval]=wavefun('db3',i);
    stairs(xval,psi)
    xlabel(i)
    pause(1)
end