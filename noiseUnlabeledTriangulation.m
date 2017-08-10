noise=0.001;
random1=rand(3,1);
random2=rand(3,1);
random3=rand(3,1);
focal_point1=6*random1/norm(random1);
focal_point2=6*random2/norm(random2);
focal_point3=6*random2/norm(random3);
A1=[eye(3) focal_point1];
A2=[eye(3) focal_point2];
A3=[eye(3) focal_point3];
twidleA1=twidleMatrix(A1);
twidleA2=twidleMatrix(A2);
twidleA3=twidleMatrix(A3);
A={twidleA1 twidleA2};
X=rand(4,1);
Y=rand(4,1);
M=X*Y'+Y*X'; M=M/M(4,4);
x1=A1*[X,Y]+normrnd(0,noise,3,2);
x2=A2*[X,Y]+normrnd(0,noise,3,2);
x3=A3*[X,Y]+normrnd(0,noise,3,2);
N1=x1(:,1)*x1(:,2)'+x1(:,2)*x1(:,1)'; vecN1 = N1(tril(true(size(N1))));
N2=x2(:,1)*x2(:,2)'+x2(:,2)*x2(:,1)'; vecN2 = N2(tril(true(size(N2))));
N3=x3(:,1)*x3(:,2)'+x3(:,2)*x3(:,1)'; vecN3 = N3(tril(true(size(N3))));
N=[vecN1 vecN2];
[kernelProjB, B]=svdTriangulation(A,N);
M0=[focal_point2;1]*[focal_point1;1]'+[focal_point1;1]*[focal_point2;1]'; M0=M0/M0(4,4);
M1=vecToSymmetricMatrix(kernelProjB(:,1));
M2=vecToSymmetricMatrix(kernelProjB(:,2));
alpha=eig(M2,M2-M1)
alpha = sym('alpha');
cubicPoly=sym2poly(det(alpha*M1+(1-alpha)*M2));
alpha=roots(cubicPoly);
beta = sym('beta');
cubicPoly=sym2poly(det(beta*M2+(1-beta)*M0));
beta=roots(cubicPoly);
mu= sym('mu');
cubicPoly=sym2poly(det(mu*M1+(1-mu)*M0));
mu=roots(cubicPoly);
nu= sym('nu');
cubicPoly=sym2poly(det(M0+nu*(M1-M2)));
nu=roots(cubicPoly);

Mres1=alpha(1)*M1+(1-alpha(1))*M2; Mres1=Mres1/Mres1(4,4);
Mres2=alpha(2)*M1+(1-alpha(2))*M2; Mres2=Mres2/Mres2(4,4);
Mres3=alpha(3)*M1+(1-alpha(3))*M2;  Mres3=Mres3/Mres3(4,4);
Mres4=alpha(4)*M1+(1-alpha(4))*M2;  Mres4=Mres4/Mres4(4,4);
M0res1=beta(1)*M2+(1-beta(1))*M0; M0res1=M0res1/M0res1(4,4);
M0res2=beta(2)*M2+(1-beta(2))*M0; M0res2=M0res2/M0res2(4,4);
M0res3=beta(3)*M2+(1-beta(3))*M0;  M0res3=M0res3/M0res3(4,4);
M0res4=beta(4)*M2+(1-beta(4))*M0;  M0res4=M0res4/M0res4(4,4);
M1res1=mu(1)*M1+(1-mu(1))*M0; M1res1=M1res1/M1res1(4,4);
M1res2=mu(2)*M1+(1-mu(2))*M0; M1res2=M1res2/M1res2(4,4);
M1res3=mu(3)*M1+(1-mu(3))*M0;  M1res3=M1res3/M1res3(4,4);
M1res4=mu(4)*M1+(1-mu(4))*M0;  M1res4=M1res4/M1res4(4,4);
M12res1=M0+nu(1)*(M1-M2); M12res1=M12res1/M12res1(4,4);
M12res2=M0+nu(2)*(M1-M2); M12res2=M12res2/M12res2(4,4);
M12res3=M0+nu(3)*(M1-M2);  M12res3=M12res3/M12res3(4,4);
M12res4=M0+nu(4)*(M1-M2);  M12res4=M12res4/M12res4(4,4);


B1=[[A1;A2] [x1(:,1);0;0;0] [0;0;0;x2(:,1)]];
B2=[[A1;A2] [x1(:,2);0;0;0] [0;0;0;x2(:,2)]];
[u,d,v]=svd(B1);
d(6,6)=0;
projB1=u*d*v';
Xrec=null(projB1);
Xrec=Xrec(1:4);
[u,d,v]=svd(B2);
d(6,6)=0;
projB2=u*d*v';
Yrec=null(projB2);
Yrec=Yrec(1:4);
Mmatched=Xrec*Yrec'+Yrec*Xrec';
Mmatched=Mmatched/Mmatched(4,4);

norm(M-Mmatched,'fro')
norm(M-Mres1,'fro')
norm(M-Mres2,'fro')
norm(M-Mres3,'fro')
norm(M-Mres4,'fro')
% norm(M-M0res1,'fro')
% norm(M-M0res2,'fro')
% norm(M-M1res1,'fro')
% norm(M-M1res2,'fro')
% norm(M-M12res1,'fro')
% norm(M-M12res2,'fro')
