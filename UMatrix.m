function U = UMatrix(H,t,n)
Id = eye(n);
h = 1.054e-34;
U = (Id - 1i*H*t/(2*h))/(Id + 1i*H*t/(2*h));