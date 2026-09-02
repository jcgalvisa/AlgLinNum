//Se va a solucionar -△u=f, con f=1, en D=[0,1]×[0,1], con u=0 en ∂D.

// Parameters
int inside = 2; 
int outside = 1;
real d=0.1;
func f = 1.0;
func g=0;       //  boundary condition function
real rate = 0; //Convergence rate

int n = 4, k = 5; // n = number of refinements, points = # of points in each borders
int nangles = 10;
real[int] Ntraingles(n), Msize(n);

//Mesh
border b1(t=0, 1){x=t; y=0; label=outside;};
border b2(t=0, 1){x=1; y=t;label=outside;};
border b3(t=1,0){x=t; y=1;label=outside;};
border b4(t=0, 1){x=0; y=t;label=outside;};

mesh Th = buildmesh(b1(k) + b2(k) + b3(k) + b4(-k));
plot(Th, wait=true);
// Loop for refinements
for(int i = 0; i < n; i++){
    // Fespace
    fespace Vh(Th, P1);
    Vh  h = hTriangle;
    Msize[i] = h[].max;
    Vh u,v;
    // Problem
    solve smoothcrack(u, v)
    = int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v))
    - int2d(Th)(f*v)
    + on(outside, u=g);
    plot(u, value=true, fill=true, wait=true);

    Ntraingles[i] = Th.nt;
    Th = trunc(Th,1,split=2); //refine mesh
    plot(Th, wait=true);
  }