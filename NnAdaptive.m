

function dx = NnAdaptive(t,x)
    gt = 1;
    kt = 1;
    s = 1;
    m = 27;
    fb = 0.26;

    pmax = 20.185;
    tb = sqrt(m)*fb/(2*pmax);

    ka1 = 1.5;
    ka2 = 2;
    ku = 2.5;
    mu = [-0.6 0 0.6]';
    phii = exp(-norm(x(1:3)-mu)^2/(2*s^2));
    phi = phii * ones(27,1);
    pt = x(4:30)'*phi;

    z1 = x(1);
    a1 = -ka1*z1;
    z2 = x(2)-a1;
    a2 = -(ka1+ka2)*z2 + z1*ka1^2;
    z3 = x(3)-a2;

    a2dot = -ka1^3*z1 + (ka1^2 + ka1*ka2 + ka2^2)*z2 -...
        (ka1 + ka2)*z3;

    zeta = 1e-6;
    db = 0.1;

    u = -ku*z3 - tanh(z3*pt/zeta)*pt - tanh(z3*a2dot/zeta)*a2dot -...
        tanh(z3*fb/(2*zeta))*fb/2-tanh(z3*db/zeta)*db;
    
    fx = x(1)^2*x(3)^2 + x(2)^4;
    gx = x(1)^2 + 2*sin(x(2))+3;
    dt = 0.1*sin(t)*sin(0.1*t);
    
    x1dot = x(2);
    x2dot = x(3);
    x3dot = fx + gx*u + dt;

    ck = 0;
     for j = 4:30
         if x(j) == 0
             ck = 1;
         end
     end

if ((norm(x(4:30)) == tb) && (z3 > 0))||((z3<0)&&(ck == 1))
         thdot = zeros(27,1);
     else
         thdot = gt*(z3*phi-kt*x(4:30));
end

dx = [x1dot;x2dot;x3dot;thdot];

