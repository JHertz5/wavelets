for i=1:10
    [phi,psi,phi_s,psi_s,xval]=wavefun('bior2.2',i);
    stairs(xval,phi_s)
    xlabel(i)
    pause(1)
end
pause(3)
for i=1:10
    [phi,psi,phi_s,psi_s,xval]=wavefun('bior2.2',i);
    stairs(xval,psi_s)
    xlabel(i)
    pause(1)
end
pause(3)
for i=1:10
    [phi,psi,phi_s,psi_s,xval]=wavefun('bior2.2',i);
    stairs(xval,phi)
    xlabel(i)
    pause(1)
end
pause(3)
for i=1:10
    [phi,psi,phi_s,psi_s,xval]=wavefun('bior2.2',i);
    stairs(xval,psi)
    xlabel(i)
    pause(1)
end