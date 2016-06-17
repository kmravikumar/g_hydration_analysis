#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "tpxio.h"

/* Places to edit in the file 

Search for "EDIT:" keywork within the file for places to edit

1) MAX_W - maximum number of water molecule sin the system

2) (X,Y,Z) simulation box size and corresponding number of cells/bins
e.g.
X_BOX_SIZE=11.0, Y_BOX_SIZE=5.5, Z_BOX_SIZE=5.0;
X_MAX_BINS=110, Y_MAX_BINS=55, Z_MAX_BINS=50;

3) Maximum Protein oxygen atoms = 500
   Maximum protein hydrogen atoms = 500

4) i - starting atom number of water oxygen-1 (C array index starts at 0)  
   Starting atom number of water oxygen

 */



const int MAX_W=11000;                                       // EDIT: Maximum Water
const float X_BOX_SIZE=11.0, Y_BOX_SIZE=5.5, Z_BOX_SIZE=5.0; // EDIT: 
const int X_MAX_BINS=110, Y_MAX_BINS=55, Z_MAX_BINS=50;      // EDIT: 

//To add two vectors
void add_vec(float a[], float b[], float c[])
  {
    c[0]= (a[0]+b[0]);
    c[1]= (a[1]+b[1]);
    c[2]= (a[2]+b[2]);
  }

//To subtract two vectors
void sub_vec(float a[], float b[], float c[])
  {
    c[0]= (a[0]-b[0]);
    c[1]= (a[1]-b[1]);
    c[2]= (a[2]-b[2]);
  }

// To get modulus of a vector
void mod_vec(float a[], float *ans)
  {
    float temp;
    temp = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);  
    *ans = sqrt(temp);
  }


// To get unit vector along a vector
void unit_vec(float a[], float unit[])
  {
    float temp; 
    temp = (a[0]*a[0])+(a[1]*a[1])+(a[2]*a[2]);  
    temp = sqrt(temp);
    unit[0] = a[0]/temp;
    unit[1] = a[1]/temp;
    unit[2] = a[2]/temp;
  }


//To find cross prod between two vectors
void cross_vec(float a[], float b[], float c[])
  {
    c[0]= (a[1]*b[2]-a[2]*b[1]);
    c[1]= (a[2]*b[0]-a[0]*b[2]);
    c[2]= (a[0]*b[1]-a[1]*b[0]);
  }


/* Put the coordinates of a[] in a 
    defined bin and return the 
    bin number as a coordinate */
void get_bin(float a[], int bin[])
  {
    int x_bin=0, y_bin=0, z_bin=0;

    x_bin = (a[0]/(X_BOX_SIZE/X_MAX_BINS));
    y_bin = (a[1]/(Y_BOX_SIZE/Y_MAX_BINS));
    z_bin = (a[2]/(Z_BOX_SIZE/Z_MAX_BINS));

    /* Due to pressure coupling the box size might 
       increase slightly > *_BOX_SIZE:
       In that case make the bin as MAX_BINS-1 */
    if(x_bin >= X_MAX_BINS){x_bin = X_MAX_BINS-1;}
    if(y_bin >= Y_MAX_BINS){y_bin = Y_MAX_BINS-1;}
    if(z_bin >= Z_MAX_BINS){z_bin = Z_MAX_BINS-1;}

    bin[0] = x_bin; 
    bin[1] = y_bin; 
    bin[2] = z_bin; 
  }


//To get dot product between two vectors
void dot_vec(float a[], float b[], float *ans)
  {
    float c[3];
    c[0]= (a[0]*b[0]);
    c[1]= (a[1]*b[1]);
    c[2]= (a[2]*b[2]);
    *ans = c[0]+c[1]+c[2];
  }

// To get the angle between vectors 
// in degrees   
void angle_bet_vec(float a[], float b[], float *theta)
  {
    float  dot_p, mod1, mod2,check;
    dot_vec(a,b,&dot_p);
    mod_vec(a,&mod1);
    mod_vec(b,&mod2);
    check = (dot_p/(mod1*mod2));
    if(check >=  1){printf("CAREFUL cos(theta) %f > 1!!",check); check =  0.999999;}
    if(check <= -1){printf("CAREFUL cos(theta) %f < 1!!",check); check = -0.999999;}
    *theta = acos(check)*180.0/3.1415;
  }

//Get H-bonding angle for ow and nei_h
float get_angle(float ow[], float nei_h[], float nei_o[])
  {
    float a[3],b[3],angle;
    a[0] = ow[0] - nei_h[0];
    a[1] = ow[1] - nei_h[1];
    a[2] = ow[2] - nei_h[2];

    b[0] = nei_h[0] - nei_o[0];
    b[1] = nei_h[1] - nei_o[1];
    b[2] = nei_h[2] - nei_o[2];

    angle_bet_vec(a, b, &angle);
    return angle;
  }  

//Get distance between points t[], b[]
float get_distance(float t[], float b[])
  {
    float x_comp, y_comp, z_comp, dist;
    float a[3];
    a[0] = t[0]; a[1] = t[1]; a[2] = t[2];

    x_comp = fabs(a[0]-b[0]);
    y_comp = fabs(a[1]-b[1]);
    z_comp = fabs(a[2]-b[2]);

    if(fabs(x_comp) >= X_BOX_SIZE/2.0 && x_comp > 0 ){x_comp = X_BOX_SIZE - x_comp;}
    if(fabs(y_comp) >= Y_BOX_SIZE/2.0 && y_comp > 0 ){y_comp = Y_BOX_SIZE - y_comp;}
    if(fabs(z_comp) >= Z_BOX_SIZE/2.0 && z_comp > 0 ){z_comp = Z_BOX_SIZE - z_comp;}

    dist = (x_comp*x_comp + y_comp*y_comp + z_comp*z_comp);
    return sqrt(dist); 
  }  

//Get sq(distance) between vectors t[], b[]
float sq_distance_vec(float t[], float b[])
  {
    float x_comp, y_comp, z_comp, dist;
    float a[3];
    a[0] = t[0]; a[1] = t[1]; a[2] = t[2];

    x_comp = (a[0]-b[0]);
    y_comp = (a[1]-b[1]);
    z_comp = (a[2]-b[2]);

    if(fabs(x_comp) >= X_BOX_SIZE/2.0 && x_comp < 0 ){a[0] = a[0] + X_BOX_SIZE;}
    if(fabs(x_comp) >= X_BOX_SIZE/2.0 && x_comp > 0 ){a[0] = a[0] - X_BOX_SIZE;}

    if(fabs(y_comp) >= Y_BOX_SIZE/2.0 && y_comp < 0 ){a[1] = a[1] + Y_BOX_SIZE;}
    if(fabs(y_comp) >= Y_BOX_SIZE/2.0 && y_comp > 0 ){a[1] = a[1] - Y_BOX_SIZE;}

    if(fabs(z_comp) >= Z_BOX_SIZE/2.0 && z_comp < 0 ){a[2] = a[2] + Z_BOX_SIZE;}
    if(fabs(z_comp) >= Z_BOX_SIZE/2.0 && z_comp > 0 ){a[2] = a[2] - Z_BOX_SIZE;}

    x_comp = fabs(a[0] - b[0]);
    y_comp = fabs(a[1] - b[1]);
    z_comp = fabs(a[2] - b[2]);
    dist = (x_comp*x_comp + y_comp*y_comp + z_comp*z_comp);
    return dist; 
  }  


// Get Orientation angle of water
// w.r.t X axis
float orientation_old(float ow[], float h1[], float h2[], float com1[], float com2[], float com3[])
  {
    float vec1[3],vec2[3];
    float angle=0.0;
    float dist1,dist2,dist3;

    sub_vec(h1,h2,vec1);
    
    dist1 = sq_distance_vec(ow,com1);
    dist2 = sq_distance_vec(ow,com2);
    dist3 = sq_distance_vec(ow,com3);

    vec2[0] = 0.0;
    vec2[1] = com1[1];
    vec2[2] = com1[2];

    if(dist2 >= dist1 && dist2 >= dist3)
      {
vec2[1] = com2[1];
vec2[2] = com2[2];
      }
    
    if(dist3 >= dist1 && dist3 >= dist2)
      {
vec2[1] = com3[1];
vec2[2] = com3[2];
      }
    
    angle_bet_vec(vec1,vec2,&angle);

    if(angle > 90)
      {
angle = 180 - angle;
      }
    return angle;
  }

// Get Orientation angle of water
// w.r.t closest collagen axes (collagen axis is parallel to x-axis) 
void orientation(float ow[], float h1[], float h2[], float com1[], float com2[], float com3[], float angle[])
  {
    float ow_com[3],vec1[3],vec2[3],vec3[3];
    float dip[3], h1h2[3], perp[3]; // These are r1,r2, and r1 X r2 vectors
    float angle1=0.0,angle2=0.0,angle3=0.0;
    float dist1,dist2,dist3;

    vec1[0] = 0.0; // collagen axis is parallel to x-axis
    vec1[1] = ow[1]-com1[1];
    vec1[2] = ow[2]-com1[2];

    vec2[0] = 0.0;
    vec2[1] = ow[1]-com2[1];
    vec2[2] = ow[2]-com2[2];

    vec3[0] = 0.0;
    vec3[1] = ow[1]-com3[1];
    vec3[2] = ow[2]-com3[2];

    mod_vec(vec1, &dist1);
    mod_vec(vec2, &dist2);
    mod_vec(vec3, &dist3);
    
    // find the closest collagen
    ow_com[0] = 0.0;
    ow_com[1] = vec1[1];
    ow_com[2] = vec1[2];

    if(dist2 <= dist1 && dist2 <= dist3)
      {
ow_com[1] = vec2[1];
ow_com[2] = vec2[2];
      }
    
    if(dist3 <= dist1 && dist3 <= dist2)
      {
ow_com[1] = vec3[1];
ow_com[2] = vec3[2];
      }

    // get the orientation of r1,r2,r3 vectors w.r.t. collagen
    sub_vec(h1,ow,vec1);
    sub_vec(h2,ow,vec2);
    add_vec(vec1,vec2,dip);    
    angle_bet_vec(dip,ow_com,&angle1);
    if(angle1 > 90)
      {
angle1 = 180 - angle1;
      }

    sub_vec(h1,h2,h1h2);    
    angle_bet_vec(h1h2,ow_com,&angle2);
    if(angle2 > 90)
      {
angle2 = 180 - angle2;
      }

    cross_vec(dip,h1h2,perp);    
    angle_bet_vec(perp,ow_com,&angle3);
    if(angle3 > 90)
      {
angle3 = 180 - angle3;
      }
    angle[0] = angle1;
    angle[1] = angle2;
    angle[2] = angle3;
  }

// Print float vector
void p_vec(float a[])
{
  printf("%f   %f   %f\n",a[0], a[1], a[2]);
}

// Print int vector
void i_vec(int a[])
{
  printf("%d   %d   %d\n",a[0], a[1], a[2]);
}


/* Main function */
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Water analysis code -Krishna(krishhere@mcbm.tamu.edu) 04/14/2010"
  };
  
  static int n=1;
  int frame=0;
  int i,j,k,a,b,c,num_atoms,r;
  float ow[MAX_W][3],h1[MAX_W][3],h2[MAX_W][3];
  float prev_ow[MAX_W][3];
  float owh1[MAX_W][3],owh2[MAX_W][3];
  float dipole[MAX_W][3],p_dipole[MAX_W][3],pp_dipole[MAX_W][3];
  float unit_dipole[MAX_W][3],unit_p_dipole[MAX_W][3], \
        unit_pp_dipole[MAX_W][3];
  float mag_temp, mag_temp1, mag_temp2, mag_temp3; 
  float prev_unit_dipole[MAX_W][3],prev_unit_p_dipole[MAX_W][3], \
        prev_unit_pp_dipole[MAX_W][3];
  float theta;
  int bin[3];  
  static int sum_density[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static int density_fluc[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  float sum_diffusion[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  float sum_dipole[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  float sum_p_dipole[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  float sum_pp_dipole[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static float sum_orientation1[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static float sum_orientation2[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static float sum_orientation3[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  float angle[3];
  static int o_occupancy[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static int prev_o_occupancy[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  static int  hydro_bond[X_MAX_BINS][Y_MAX_BINS][Z_MAX_BINS];
  int l,m,p,e,f,g;
  int min_x,max_x,min_y,max_y,min_z,max_z;
  int water_num, h;
  float nei_h1[3],nei_h2[3];;
  float h_dist1,h_dist2, h_angle1, h_angle2;
  float vec1[3],vec2[3],vec3[3];
  float com1[3],com2[3],com3[3];
  float x_sum=0,y_sum=0,z_sum=0;
  int num_of_O=0, prot_O[500]; /* EDIT: Maximum of 500 Protein Oxygen atoms */
  int num_of_H=0, prot_H[500]; /* EDIT: Maximum of 500 Protein Hydrogen atoms */
  static float prot_O_coord[500][3], prot_H_coord[500][3]; /* EDIT: Maximum of 500 Protein Oxygen/Hydrogen atoms */

  /* initialize center of mass of 3 collagen peptides */
  for(i=0;i<3;i++)
    {
      com1[i]=0;
      com2[i]=0;
      com3[i]=0;      
    }

  /* initialize arrays to 0 */
  for(i=0;i<X_MAX_BINS;i++)
    {
      for(j=0;j<Y_MAX_BINS;j++)
{
  for(k=0;k<Z_MAX_BINS;k++)
    {
      density_fluc[i][j][k] = 0;
      sum_density[i][j][k] = 0;
      sum_diffusion[i][j][k] = 0.0;
      sum_dipole[i][j][k] = 0.0;
      sum_p_dipole[i][j][k] = 0.0;	
      sum_pp_dipole[i][j][k] = 0.0;
      sum_orientation1[i][j][k] = 0.0;
      sum_orientation2[i][j][k] = 0.0;
      sum_orientation3[i][j][k] = 0.0;
      o_occupancy[i][j][k] = -1;
      prev_o_occupancy[i][j][k] = -1;
      hydro_bond[i][j][k] = 0;
    }
}
    }

  /* initialize Prot_*_arrays to 0 */
  for(j=0;j<500;j++)
    {
      for(i=0;i<3;i++)
{
  prot_O_coord[j][i]=0.0;
  prot_O_coord[j][i]=0.0;
  prot_O_coord[j][i]=0.0;      
  prot_H_coord[j][i]=0.0;
  prot_H_coord[j][i]=0.0;
  prot_H_coord[j][i]=0.0;     
}
    }

  /* Extra arguments - but note how you always get the begin/end
• options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  
  t_topology top;
  int        ePBC;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status;
  int        flags = TRX_READ_X;

  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
• calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* We don't need any topology information to write the coordinates,
• but to show how it works we start by writing the name and
• charge of the selected atom. It returns a boolean telling us
• whether the topology was found and could be read
   */
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
  sfree(xtop);

  /*printf("Atom name: %s\n",*(top.atoms.atomname[n]));
    printf("Atom charge: %f\n",top.atoms.atom[n].q); */
  
  /* The first time we read data is a little special */
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

  /* frame number */
  frame = 0;
   
  /* This is the main loop over frames */
    do {
    /* coordinates are available in the vector fr.x
• you can find this and all other structures in
`    * the types directory under the gromacs include dir.
• Note how flags determines wheter to read x/v/f!
     */

      /* initialize cell occupancy arrays to 0 for each iteration */
      for(i=0;i<X_MAX_BINS;i++)
{
  for(j=0;j<Y_MAX_BINS;j++)
    {
      for(k=0;k<Z_MAX_BINS;k++)
{
  // If previous occupancy of cell is different from current
  // occupancy increase density fluctuation by 1
  if(prev_o_occupancy[i][j][k] != o_occupancy[i][j][k]) 
    {
      density_fluc[i][j][k] += 1;
    }
  prev_o_occupancy[i][j][k] = o_occupancy[i][j][k];  
  o_occupancy[i][j][k] = -1;
}
    }
}

    /* Calculate COM of the 3 collagen peptides */    
    if(frame==0)
      {
num_atoms = 0;
x_sum=0; y_sum=0; z_sum=0;
for(i=0;i<927;i++) //EDIT: 1st molecule index
  {
    num_atoms++;
    x_sum += fr.x[i][XX];
    y_sum += fr.x[i][YY];
    z_sum += fr.x[i][ZZ];
  }
com1[0] = x_sum/num_atoms;
com1[1] = y_sum/num_atoms;
com1[2] = z_sum/num_atoms;
	
num_atoms = 0;
x_sum=0; y_sum=0; z_sum=0;
for(i=927;i<1875;i++) //EDIT: 2nd molecule index
  {
    num_atoms++;
    x_sum += fr.x[i][XX];
    y_sum += fr.x[i][YY];
    z_sum += fr.x[i][ZZ];
  }
com2[0] = x_sum/num_atoms;
com2[1] = y_sum/num_atoms;
com2[2] = z_sum/num_atoms;
	
num_atoms = 0;
x_sum=0; y_sum=0; z_sum=0;
for(i=1875;i<2822;i++) //EDIT: 3rd molecule index
  {
    num_atoms++;
    x_sum += fr.x[i][XX];
    y_sum += fr.x[i][YY];
    z_sum += fr.x[i][ZZ];
  }
com3[0] = x_sum/num_atoms;
com3[1] = y_sum/num_atoms;
com3[2] = z_sum/num_atoms;

printf("\n COMs of 3 peptides = ");
p_vec(com1); p_vec(com2); p_vec(com3);

    /* Store Protein Oxygen,Hydrogen coordinates for H-bond distance checking */
// Number of oxygens and hydrogens 
// Initialize to 0
num_of_O = 0;
num_of_H = 0;
for(i=0; i<2823; i++) //EDIT: number of protein atoms
  {
    if(strcmp(*top.atoms.atomname[i],"O")==0  \
    || strcmp(*top.atoms.atomname[i],"OG2")==0 )
      {
prot_O[num_of_O] = i;
num_of_O += 1;
      }
    if(strcmp(*top.atoms.atomname[i],"HN")==0  \
    || strcmp(*top.atoms.atomname[i],"HG2")==0 )
      {
prot_H[num_of_H] = i;
num_of_H += 1;
      }
  }
printf("\n number of Protein Oxygens = %d \n", num_of_O);
printf("number of Protein Hydrogens = %d \n", num_of_H);
      }
	    
    /* assign bin numbers to Proteins Oxygen and Hydrogen atoms */
    if(frame>=0)
      {
// prot_O_coord[r] has coordinates of protein O atom r
for(r=0;r<num_of_O;r++)
  {
    prot_O_coord[r][0] = fr.x[prot_O[r]][XX];
    prot_O_coord[r][1] = fr.x[prot_O[r]][YY];
    prot_O_coord[r][2] = fr.x[prot_O[r]][ZZ];
  }
for(r=0;r<num_of_H;r++)
  {
    prot_H_coord[r][0] = fr.x[prot_H[r]][XX];
    prot_H_coord[r][1] = fr.x[prot_H[r]][YY];
    prot_H_coord[r][2] = fr.x[prot_H[r]][ZZ];
  }
      }

    /* For all water oxygen atoms */
    j = 0; // j is the water oxygen number

    // EDIT: i - starting atom number of water oxygen-1 (C array index starts at 0)
    // increment by 3 for 3 water atoms (ow,h1,h2)
    for(i=2823;i<(top.atoms.nr);i=i+3) 
      {
// assign coordinates of j th oxygen
// and hydrogen atoms to ow[j], h1[j], h2[j]
ow[j][0] = fr.x[i][XX];   
ow[j][1] = fr.x[i][YY];   
ow[j][2] = fr.x[i][ZZ]; 
  
h1[j][0] = fr.x[i+1][XX]; 
h1[j][1] = fr.x[i+1][YY]; 
h1[j][2] = fr.x[i+1][ZZ];   

h2[j][0] = fr.x[i+2][XX]; 
h2[j][1] = fr.x[i+2][YY]; 
h2[j][2] = fr.x[i+2][ZZ];   


// Calculate unit vectors along 
// 1. (r1) dipole
// 2. (r2) perpendicular to dipole (p_dipole)
// 3. (r3) perpendicular to both the above (pp_dipole)
sub_vec(h1[j],ow[j],owh1[j]);
sub_vec(h2[j],ow[j],owh2[j]);
add_vec(owh1[j],owh2[j],dipole[j]);
sub_vec(owh1[j],owh2[j],p_dipole[j]);

unit_vec(dipole[j],unit_dipole[j]);
unit_vec(p_dipole[j],unit_p_dipole[j]);
cross_vec(unit_dipole[j],unit_p_dipole[j],unit_pp_dipole[j]);

// get water oxygen bin number
get_bin(ow[j],bin);
a = bin[0];
b = bin[1];
c = bin[2];

// occupancy of bin a,b,c is the water number (j) or -1
o_occupancy[a][b][c] = j;

// increase density count of bin a,b,c by 1
sum_density[a][b][c] += 1;

// find orientation angle of r1, r2,r3 vectors 
// w.r.t 3 collagen com
orientation(ow[j],h1[j],h2[j],com1,com2,com3,angle);

// Sum the orientation angles in each bin
sum_orientation1[a][b][c] += angle[0];
sum_orientation2[a][b][c] += angle[1];
sum_orientation3[a][b][c] += angle[2];

/* From the second frame onwards, compare r1,r2,r3 vectors 
// for each water from the two consecutive frames and find their
// rotation angles.
// Sum rotation angles of each vector to the previous bin
// containing that water oxygen atom. */
if(frame>=1)
  {
    mag_temp = sq_distance_vec(ow[j],prev_ow[j]);

    dot_vec(unit_dipole[j],prev_unit_dipole[j], &mag_temp1);
    if(mag_temp1 >  1.00){mag_temp1 =  0.999999;}
    if(mag_temp1 < -1.00){mag_temp1 = -0.999999;}
    mag_temp1 = acos(mag_temp1);

    dot_vec(unit_p_dipole[j],prev_unit_p_dipole[j], &mag_temp2);
    if(mag_temp2 >  1.00){mag_temp2 =  0.999999;}
    if(mag_temp2 < -1.00){mag_temp2 = -0.999999;}
    mag_temp2 = acos(mag_temp2);

    dot_vec(unit_pp_dipole[j],prev_unit_pp_dipole[j], &mag_temp3);
    if(mag_temp3 >  1.00){mag_temp3 =  0.999999;}
    if(mag_temp3 < -1.00){mag_temp3 = -0.999999;}
    mag_temp3 = acos(mag_temp3);
	    
    sum_diffusion[a][b][c] += mag_temp;
    sum_dipole[a][b][c]    += mag_temp1*mag_temp1;
    sum_p_dipole[a][b][c]  += mag_temp2*mag_temp2;
    sum_pp_dipole[a][b][c] += mag_temp3*mag_temp3;

  }
// increment water oxygen
j++;
      }

    
    /* Again for all water atoms 
       Second loop for each water oxygen atom 
       Variable j has the total number of oxygen atoms */
    // Hydogen bonds are calculated in the current frame

    for(k=0;k<j;k=k++)
      {

// prev_* arrays store data/coordinates from the 
// previous frame
prev_ow[k][0] = ow[k][0]; 
prev_ow[k][1] = ow[k][1]; 
prev_ow[k][2] = ow[k][2];
prev_unit_dipole[k][0] = unit_dipole[k][0];
prev_unit_dipole[k][1] = unit_dipole[k][1];
prev_unit_dipole[k][2] = unit_dipole[k][2];
prev_unit_p_dipole[k][0] = unit_p_dipole[k][0];
prev_unit_p_dipole[k][1] = unit_p_dipole[k][1];
prev_unit_p_dipole[k][2] = unit_p_dipole[k][2];
prev_unit_pp_dipole[k][0] = unit_pp_dipole[k][0];
prev_unit_pp_dipole[k][1] = unit_pp_dipole[k][1];
prev_unit_pp_dipole[k][2] = unit_pp_dipole[k][2];

// a,b,c has bin numbers of water oxygen k
get_bin(ow[k],bin);
a = bin[0];
b = bin[1];
c = bin[2];

// search for next oxygen atom 
// in a cube surrounding the current cell
// NOTE: search for water oxygens not hydrogens
// So the the cube has to be bigger (+5,-5)
min_x = a-5; max_x = a+5;
min_y = b-5; max_y = b+5;
min_z = c-5; max_z = c+5;

/* Count hydrogen bonding with Protein atoms*/
// For each protein oxygen atom fing distance between 
// hydrogens of water k.
// If distance <0.24 then increment H-bonding of 
// bin a,b,c by 1
for(l=0;l<num_of_O;l++)
  {
    h_dist1 = get_distance(h1[k],prot_O_coord[l]);
    h_dist2 = get_distance(h2[k],prot_O_coord[l]);
    if(h_dist1 < 0.24 || h_dist2 < 0.24)
      {
hydro_bond[a][b][c] += 1;
      }
  }

// Do the same as above for protein hydrogen atoms
for(l=0;l<num_of_H;l++)
  {
    h_dist1 = get_distance(ow[k],prot_H_coord[l]);
    if(h_dist1 < 0.24 )
      {
hydro_bond[a][b][c] += 1;
      }
  }

// for each cell loop through the neighbouring cells
// and check for hydrogen bonds.
// If the cells are at the edges, then use 
// periodic boundary condition to find the neighbours
// NOTE: Using periodic boundary condition is not accurate
// when box size fluctuates under pressure coupling
for(l=min_x;l<=max_x;l++)
  {
    e = l;
    if(l < 0){ e = e + X_MAX_BINS; } 
    if(l >= X_MAX_BINS){ e = e - X_MAX_BINS; } 
    for(m=min_y;m<=max_y;m++)
      {
f = m;
if(m < 0){f = f + Y_MAX_BINS;}
if(m >= Y_MAX_BINS){f = f - Y_MAX_BINS;}
for(p=min_z;p<=max_z;p++)
  {
    g = p;
    if(p < 0){g = g + Z_MAX_BINS;}
    if(p >= Z_MAX_BINS){g = g - Z_MAX_BINS;}
    if(o_occupancy[e][f][g] != -1 )
      {
water_num = o_occupancy[e][f][g];
nei_h1[0] = h1[water_num][0];
nei_h1[1] = h1[water_num][1];
nei_h1[2] = h1[water_num][2];
nei_h2[0] = h2[water_num][0];
nei_h2[1] = h2[water_num][1];
nei_h2[2] = h2[water_num][2];
h_dist1 = get_distance(ow[k],nei_h1);
h_dist2 = get_distance(ow[k],nei_h2);
//h_angle1 = get_angle(ow[k],nei_h1,nei_o);
//h_angle2 = get_angle(ow[k],nei_h2,nei_o);
h_angle1=0.0;
h_angle2=0.0;

if(h_dist1 <= 0.24 && h_angle1 <=40.0 && k!= water_num )
  {
    hydro_bond[a][b][c] += 1;
    hydro_bond[e][f][g] += 1;
  }
if(h_dist2 <= 0.24 && h_angle2 <=40.0 && k!= water_num )
  {
    hydro_bond[a][b][c] += 1;
    hydro_bond[e][f][g] += 1;
  }
      }
  }
      }
  }
	    
      }

    frame += 1;
    //if(fr.time>=2000.0){break;}
  } while(read_next_frame(status,&fr));   /* This is the main loop over frames */


  /* Print total water in the system */
  printf("\n Total Water = %d \n",j);


  /* Write Data to File */
  FILE *fpw1  ;
  fpw1 = fopen("final_data.dat","w");
  
  for(i=0;i<X_MAX_BINS;i++)
    {
      for(j=0;j<Y_MAX_BINS;j++)
{
  for(k=0;k<Z_MAX_BINS;k++)
    {
      fprintf(fpw1,"%4d %4d %4d %5d %11.5f %11.5f %11.5f %11.5f %13.5f %13.5f %13.5f %6d %6d \n", \
      i,j,k,sum_density[i][j][k],sum_diffusion[i][j][k], \
      sum_dipole[i][j][k], sum_p_dipole[i][j][k], sum_pp_dipole[i][j][k], \
      sum_orientation1[i][j][k], sum_orientation2[i][j][k], sum_orientation3[i][j][k], \
      hydro_bond[i][j][k], density_fluc[i][j][k]);
    }
}
    }
  
  fclose(fpw1);
  return 0;
} 

