//********** 3-dimensional registration of particles derived from two images ************************************
//
// Written by Hiroshi Koyama in National Institute for Basic Biology, Japan, 2022/02/03
// Revised by Hiroshi Koyama, 2022/07/21
// hkoyama@nibb.ac.jp
//
// <Description>
// Read a text file (input_xyz_registration.csv) containing xyz-coordinates of multiple particles,
// 	where the particles are derived from two sample images (e.g. one before fixation, the other after fixation).
// Here, the two sample images called image-before and image-after.
// Among the particles, several particles are landmark particles, and the other particles are not landmark particles.
// Users should set so that each landmark particle has the same ID between the before- and after-images in the text file,
// 	resulting in multiple sets of paired landmark particles.
// The particles other than the landmark particles should have different IDs between the two images.
// 3-dimesional registration is performed toward the xyz-coordinates of the landmark particles,
// 	where the summation of distances between the paired landmark particles between the before- and after-images is minimized
// 	by a Monte Carlo algorithm (MCMC: Markov chain Monte Carlo method).
// 	**** Relatd result files = 	output_rotation_results.txt
// 								output_bad_results_rotation.txt

// If the registration is successful, the xyz-coordinates of the landmark particles after the registration are saved
// 	with the distances between the paired landmark particles. 
// 	**** Related result file = 	output_registered_particles_landmark.txt 

// Then, the particles other than the landmark particles are focused on.
// For each particle in the before-image, one or three nearest neighbor particles are searched for the after-image, 
// 	which is/are the candidate(s) of the pair of the particle in the before-image.
// These candidate particles are saved with their xyz-coordinates and with the distances from the particle in the before-image.
// A similar process is also performed for each particle in the after-image toward the before-image. 
// 	**** Result file = 	output_registered_particles_1neighbors.txt
// 						output_registered_particles_3neighbors.txt
// 						output_registered_particles_all.txt

// For trouble shooting.
// The minization process is based on MCMC where kBT and end_iteration parameters affect the success/failure of minimization.
// If minimization results are not satisfactory, the above two parameters may improve the results. 
// Also, see the section of "Setting of initial angular configuration number" below. 

// For 3D visualization.
// The result of the 3D registration can be visualized through Macro_particle_drawing_220721_02.ijm and/or Macro_3D_image_rotation_220722_2.ijm.
// To achieve good visualization, multiple parameter setting about scales of images will be required.
// To avoid the complicated parameter setting and the troubles about the scales,
// I recommend that you set the same scale among image-before, image-after, and input_xyz_registration.csv.
//		The scales include not only xy scales but also z scales, i.e. um/pixel and um/voxel. 

// A tip of obtained angles
// 3-dimensional rotation is defined by three angles, whereas 2-dimensional rotation is defined by one angle.
// For 3-dimensional rotation, different sets of three angles can achieve the same rotation, 
//  which is contrast to 2-dimensional situations.
// Therefore, we may obtain different sets of three angles from different initial sets of the three angles, 
//  all of which correctly register the particles.

// <Requirement>
// Fiji, version 1.53 or maybe later.
// The name of the text file is input_xyz_registration.csv with a given format.
// The number of landmark particles should be larger than 3, and the registration may work well in the case that the number is ~10.
// We experinced systems with the number of all particles in a before-image = 50-100, the number of landmark particles = 8~10.

// <Running time>
// MacBook air (CPU: Core i7, Dual core, 1.7GHz; Memory: 8 GB 1600 MHz DDR3)
// In the case that number_angle = 3 and landmark pariticles = 9, the running time was ~50sec.

// <Reference>
// Koyama et al. "An ImageJ-based tool for three-dimensional registration between different types of microscopic images"
// Cite the reference above when you use this macro for any publication.
//
//*****************************************************************************************************************/

//********** Parameters to be defined by users ********************************************************************

print("Please set several values in 'Parameters to be defined by users' on the macro code.\n");

//** Scaling of xyz-coordinates: if the values = 1.0, the xyz-coordinates from the input file are not modified. 
//In spite of these parameters, I recommend that scaling of xyz-coordinates is performed toward the original image but not toward the text file.
//This is because the scaling toward the text file cause troubles for visualization of the results using 
// "Macro_particle_drawing_220721_02.ijm and "Macro_3D_image_rotation_220722_2.ijm": incorrect visualizetion may occur. 
//for before-image
x_scale_before = 1.0;
y_scale_before = 1.0;
z_scale_before = 1.0;
//for after-image
x_scale_after = 1.0;
y_scale_after = 1.0;
z_scale_after = 1.0;

//** Setting a threshold value for judging successful minimization in MCMC
// This value is given as mean distance between paired particles in the x-scale of before-image.
// e.g. if x-length of a system is 100.0 with 10 linearly aligned particles, the threshold value = ~10.0 may be good.
distance_threshold = 10.0;

//** Setting of initial angular configuration number
// Segment number of 3-angles for initial values of angle search.
// First, try number_angle = 1 or 3.
// Then, if you cannot obtain a good result, increase number_angle; e.g. 5, 7, 9. (Odd numbers are preferable.)
// If number_angle = 3, 3^3=27 patterns of the intial angles are applied. It may takes a few minutes for finishing one run.
// If number_angle = 9, 9^3=729 patterns, thus, time consuming. It may takes ~1 hour for finishing one run.
// If you cannot find good results under number_angle=9, maybe something wrong.
// For trouble shooting, see the section "For trouble shooting" in <Description> above.
number_angle = 3;

//*******************************************************************************************************************


//***** input file name *****
dir_input = getDirectory("Choose a Directory ");
input_file = dir_input+"input_xyz_registration.csv";

//***** Setting of random number *****
random("seed", 100);

//***** global variables ***********************************************
var c_num=1000;
var kBT=1.0000;				//In our samples, 1.0~3.0 is good.		
var mag=1.0;				//In our samples, 0.01~1.0 is good.
var end_iteration=30000;	//In our samples, 30000 is sufficient, but 10000 is not good.

var rfx, rfy, rfz;
var max_x, max_y, max_z, min_x, min_y, min_z;
var buf_x, buf_y, buf_z;
var x_scaling, y_scaling, z_scaling;
var x_scaling2, y_scaling2, z_scaling2;
var cent_displace_x, cent_displace_y, cent_displace_z;
var rot_displace_phi, rot_displace_theta, rot_displace_psi;
var phi, psi, theta;
var centroid_x1, centroid_y1, centroid_z1, centroid_x2, centroid_y2, centroid_z2;
var Xm, Ym, Zm, Xn, Yn, Zn, Xo, Yo, Zo, Xq, Yq, Zq;
var Xp, Yp, Zp, c_pro_cost1, c_pro_cost2, min_c_pro_cost, initial_cost, c1, c2, c3, s1, s2, s3;
var random1, frequency, error_D;
var min_dist2, dist2, min_dist, dist;
var xyz_err1, mag_xyz;
var longest_axis, scaling_val;
var dist_threshold, dist_thres_1;

var cell_order, cell_num_t;
var cell_share1, cell_pair, paired_id;
var max_cell_number, cell_number;
var buf_t, buf_t2, buf_id;
var end_t;
var t, tt, g, i, ii, j, p1, p2, p3, count0, count1, count2, count3, count4, count5, count7, count8;
var ang, min_iter;
var min_id, id_chofuku, id_chofk_num;
var min_id_order;
var d_ang;
var landmark_num;

rfx = newArray(2*c_num);
rfy = newArray(2*c_num);
rfz = newArray(2*c_num);
Xm = newArray(c_num);
Ym = newArray(c_num);
Zm = newArray(c_num);
Xn = newArray(c_num);
Yn = newArray(c_num);
Zn = newArray(c_num);
Xo = newArray(c_num);
Yo = newArray(c_num);
Zo = newArray(c_num);
Xq = newArray(c_num);
Yq = newArray(c_num);
Zq = newArray(c_num);
rot_displace_phi = newArray(2);
rot_displace_theta = newArray(2);
rot_displace_psi = newArray(2);
phi = newArray(3);
psi = newArray(3);
theta = newArray(3);
min_dist2 = newArray(3);

cell_order = newArray(2*c_num);
cell_num_t = newArray(c_num);
cell_share1 = newArray(2*c_num);
cell_pair = newArray(c_num);
paired_id = newArray(c_num);
id_chofuku = newArray(c_num);
id_chofk_num = newArray(c_num);
min_id = newArray(3);

//*************************************************************************

//***** basic parameter setting
x_scaling = x_scale_before;
y_scaling = y_scale_before;
z_scaling = z_scale_before;
x_scaling2 = x_scale_after;
y_scaling2 = y_scale_after;
z_scaling2 = z_scale_after;
d_ang = number_angle;
dist_thres_1 = distance_threshold;

//***** input file reading *****
input_file_reading (input_file);

//***** xyz scale revision *****/
xyz_scaling ();

//***** extraction of landmark particles (t=0) and (t=1) *****/
landmark_particle_extraction ();

//***** modification of xyz: centroid-conserved = translational motion removed *****/
translational_motion_removal ();

//***** modification of xyz: centroids are moved to be (0,0,0). *****/
translation_movement_to_origin ();

//***** Setting of manitude of xyz-error related to MCMC ******/
// This value affects the reject/accept process in MCMC.
// The value is automatically detemined by the system size and particle numbers.
// Basically, the longest axis length in the system is fixed to be 100.0 in the xyz_scaling () function.
// Then, if the particle number is ~10, mag_xyz=0.03 may work well.
// In MCMC, the value of the xyz-error affects the success/failure of the minimization.
// In systems which can exhibit larger xyz-error values, mag_xyz should be set to be smaller.
// By contrast, systems with smaller xyz-error values, mag_xyz should be set to be larger.
// If you want to set the mag_xyz value by your hand, you can do that.
mag_xyz=0.03;
mag_xyz=mag_xyz*(100.0/(longest_axis*scaling_val))*(10.0/parseFloat(landmark_num));


//***** MCMC for xyz_error minimization by rotating after-image particles  *****/
//* output files */
output_1 = dir_input+"output_rotation_results.txt";
output_2 = dir_input+"output_bad_results_rotation.txt";

File.saveString("List of successful MCMC", output_1);
File.append("	xyz-error: improvement of xyz-error by MCMC", output_1);
File.append("	mean-distance: mean distance between paired landmark particles after xyz-error minimization", output_1);
File.append("		1st mean-distance = distance in the original scale corresponding to x-scale in before-image", output_1);
File.append("		2nd mean-distance = distance in the revised scale used for MCMC", output_1);
File.append("	angle: three angles for rotation", output_1);

File.saveString("List of failed MCMC", output_2);
File.append("	xyz-error: improvement of xyz-error by MCMC", output_2);
File.append("	mean-distance: mean distance between paired landmark particles after xyz-error minimization", output_2);
File.append("		1st mean-distance = distance in the original scale corresponding to x-scale in before-image", output_2);
File.append("		2nd mean-distance = distance in the revised scale used for MCMC", output_2);
File.append("	angle: three angles for rotation", output_2);


//* Multiple initial rotation angles are generated according to d_ang parameter.
// For each initial rotation angles, MCMC-based minimization of xyz error is performed.
// Among them, the best result (i.e. xyz error is smallest) is selected.
min_dist = 1000000000.0;
t=1;	//t=1 means the second (after-) image
for (p1=0; p1<d_ang; p1++){
	phi[0] = 2.0*PI*parseFloat(p1)/parseFloat(d_ang);

	for (p2=0; p2<d_ang; p2++){
		phi[1] = 2.0*PI*parseFloat(p2)/parseFloat(d_ang);

		for (p3=0; p3<d_ang; p3++){
			phi[2] = 2.0*PI*parseFloat(p3)/parseFloat(d_ang);

			//MCMC for xyz minimization
			xyz_minimization ();

			//Select good conditions which give small xyz error:
			//e.g. mean distance between paired particles is less than a threshold.
			distance_error_evaluation ();
			dist2 = dist/(x_scaling*scaling_val);
			if (dist<dist_threshold){
				str3 = toString(t)+": "+toString(min_iter)+" xyz-error("+toString(initial_cost)+" -> "+toString(min_c_pro_cost)+") mean-distance("+toString(dist2)+", "+toString(dist)+") angle("+toString(theta[0])+", "+toString(theta[1])+", "+toString(theta[2])+")";
				File.append(str3, output_1);

				//The angles for the best condition. i.e. a condition giving the smallest xyz error
				if (dist<min_dist){
					min_dist = dist;
					rot_displace_phi[t]=theta[0];
				    rot_displace_theta[t]=theta[1];
				    rot_displace_psi[t]=theta[2];
				}
			}
			else{
				str3 = toString(t)+": "+toString(min_iter)+" xyz-error("+toString(initial_cost)+" -> "+toString(min_c_pro_cost)+") mean-distance("+toString(dist2)+", "+toString(dist)+") angle("+toString(theta[0])+", "+toString(theta[1])+", "+toString(theta[2])+")";
				File.append(str3, output_2);
			}
		}
	}
}

//The best rotatitional angles are saved.
str3 = "best rotational angles: "+toString(rot_displace_phi[t])+", "+toString(rot_displace_theta[t])+", "+toString(rot_displace_psi[t]);
File.append(str3, output_1);

//***** Rotational angle modification by using the best rotational matrix calculated above *****/
rotation_angle_modification ();	

//**** xyz scale revision ****/
//In xyz_scaing () function, some modifications were performed for xyz-coordinates for the purpose of the minization process.
//To fairly compare the output xyz-coordinates with the distance_threshold parameter defined by users, the xyz-scales are revised.
xyz_scale_reverse ();

//***** output of pairing between before- and after-images *****/
//1-neighbors or 3-neighbors are saved.
//Two types of output: before-image as a standard or afterimage as a standard
output_pairing_particles ();

//***** output of xyz coordinates of landmark or all particles after 3D registration, and their xyz error(=distance)  *****/
output_xyz_coordinate_after_rotation ();

//***** output of condition *****/
output_conditions (input_file);

//******************* end of macro ********************************************


//******************* functions ***********************************************

function rotation_angle_modification ()
{
	t=1;
	for (i=0; i<=cell_num_t[t]; i++){
		Xm[i]=rfx[t*c_num+i];
	    Ym[i]=rfy[t*c_num+i];
	    Zm[i]=rfz[t*c_num+i];
	}

	tt=1;
	//* calculation of a rotational center = centroid */
	count0=0;
	centroid_x2=0.0;
    centroid_y2=0.0;
    centroid_z2=0.0;
	for (i=0; i<=cell_num_t[tt]; i++){
		if (cell_share1[tt*c_num+i]==1){
			count0++;
			centroid_x2=centroid_x2+rfx[tt*c_num+i];
			centroid_y2=centroid_y2+rfy[tt*c_num+i];
			centroid_z2=centroid_z2+rfz[tt*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x2=centroid_x2/parseFloat(count0);
		centroid_y2=centroid_y2/parseFloat(count0);
		centroid_z2=centroid_z2/parseFloat(count0);
	}

	//* xyz-coordinates transformed to rotational center */
	for (i=0; i<=cell_num_t[t]; i++){
		Xn[i]=Xm[i]-centroid_x2;
	    Yn[i]=Ym[i]-centroid_y2;
	    Zn[i]=Zm[i]-centroid_z2;
	}

	//* rotation angle */
	psi[0]=rot_displace_phi[tt];
    psi[1]=rot_displace_theta[tt];
    psi[2]=rot_displace_psi[tt];
	c1=cos(psi[0]);
    c2=cos(psi[1]);
    c3=cos(psi[2]);
	s1=sin(psi[0]);
    s2=sin(psi[1]);
    s3=sin(psi[2]);

	//* rotation: using rotation matrix with yaw-pitch-roll (not Euler angles) */
	for (i=0; i<=cell_num_t[t]; i++){
		Xo[i]=(Xn[i]*c1*c2)+(Yn[i]*(c1*s2*s3-s1*c3))+(Zn[i]*(c1*s2*c3+s1*s3));
		Yo[i]=(Xn[i]*s1*c2)+(Yn[i]*(s1*s2*s3+c1*c3))+(Zn[i]*(s1*s2*c3-c1*s3));
		Zo[i]=(Xn[i]*(-s2))+(Yn[i]*c2*s3)+(Zn[i]*c2*c3);
	}

	//* xyz-coordinates transformed to the original xyz-space rather than rotational center */
	for (i=0; i<=cell_num_t[t]; i++){
		Xm[i]=Xo[i]+centroid_x2;
	    Ym[i]=Yo[i]+centroid_y2;
	    Zm[i]=Zo[i]+centroid_z2;
	}

	for (i=0; i<=cell_num_t[t]; i++){
		rfx[t*c_num+i]=Xm[i];
	    rfy[t*c_num+i]=Ym[i];
	    rfz[t*c_num+i]=Zm[i];
	}
}

function xyz_minimization ()
{
	psi[0]=phi[0];
    psi[1]=phi[1];
    psi[2]=phi[2];
	theta[0]=phi[0];
    theta[1]=phi[1];
    theta[2]=phi[2];
	xyz_err1 = 0.0;

	//* initial calculation of xyz error */
	for (i=0; i<=count5; i++){
		Xp=Xn[i]-Xm[i];
	    Yp=Yn[i]-Ym[i];
	    Zp=Zn[i]-Zm[i];
		Xq[i]=Xn[i];
	    Yq[i]=Yn[i];
	    Zq[i]=Zn[i];
		xyz_err1 = xyz_err1 + Xp*Xp + Yp*Yp + Zp*Zp;
	}
	c_pro_cost1 = xyz_err1;
	min_c_pro_cost=c_pro_cost1;
	initial_cost=c_pro_cost1;
	min_iter=0;

	//* MCMC for phi to minimize c_pro_cost */
	if (count5>0){	//count5 = number of landmark particles
		for (g=0; g<=end_iteration; g++){

			random1=random;
			ang=parseInt(random1*3.0-0.5);	//In ImageJ macro, parseInt(2.4) -> 2, parseInt(2.5) -> 3.
			if (ang<=-1){	ang=0;	}
			if (ang>=3){	ang=2;	}	//necessary
			random1=random;
			random1=mag*(random1-0.5);
			psi[ang]=psi[ang]+random1;
			if (psi[ang]>PI*2.0){
				psi[ang] = psi[ang] - PI*2.0;
			}
			if (psi[ang]<0.0){
				psi[ang] = psi[ang] + PI*2.0;
			}
			c1=cos(psi[0]);
		    c2=cos(psi[1]);
		    c3=cos(psi[2]);
			s1=sin(psi[0]);
		    s2=sin(psi[1]);
		    s3=sin(psi[2]);

			//* rotation: using rotation matrix with yaw-pitch-roll (not Euler angles) */
			for (i=0; i<=count5; i++){	//count5 = number of landmark particles
				Xo[i]=(Xn[i]*c1*c2)+(Yn[i]*(c1*s2*s3-s1*c3))+(Zn[i]*(c1*s2*c3+s1*s3));
				Yo[i]=(Xn[i]*s1*c2)+(Yn[i]*(s1*s2*s3+c1*c3))+(Zn[i]*(s1*s2*c3-c1*s3));
				Zo[i]=(Xn[i]*(-s2))+(Yn[i]*c2*s3)+(Zn[i]*c2*c3);
			}

			//* MCMC: minimize xyz_error */
			xyz_err1=0.0;
			for (i=0; i<=count5; i++){
				Xp=Xo[i]-Xm[i];
			    Yp=Yo[i]-Ym[i];
			    Zp=Zo[i]-Zm[i];
				xyz_err1 = xyz_err1 + Xp*Xp + Yp*Yp + Zp*Zp;	
			}
			c_pro_cost2 = xyz_err1;	
		
			error_D = c_pro_cost2 - c_pro_cost1;
			if (error_D<=0.0){
				frequency=1.0;
			}
			else{
				frequency=exp(-error_D*mag_xyz/kBT);
			}
			random1=random;
			if (random1>frequency){	psi[ang]=phi[ang];	}	//rejected
			else{											//accepted
				phi[ang]=psi[ang];
				c_pro_cost1=c_pro_cost2;
				if (c_pro_cost1<=min_c_pro_cost){
					min_c_pro_cost=c_pro_cost1;
					theta[0]=phi[0]; 
					theta[1]=phi[1];
				    theta[2]=phi[2];
					min_iter=g;
					for (i=0; i<=count5; i++){
						Xq[i]=Xo[i];
					    Yq[i]=Yo[i];
					    Zq[i]=Zo[i];
					}
				}
			}	
		}
	}
	else{	//If there is no landmark particles defined, do nothing.
		theta[0]=0.0;
	    theta[1]=0.0;
	    theta[2]=0.0;
    }
}

function distance_error_evaluation ()
{
	dist=10000000.0;
	if (count5>0){
		dist=0;		

		c1=cos(theta[0]); 
		c2=cos(theta[1]); 
		c3=cos(theta[2]);
		s1=sin(theta[0]); 
		s2=sin(theta[1]);
	    s3=sin(theta[2]);

		//* rotation: using rotation matrix with yaw-pitch-roll (not Euler angles) */
		for (i=0; i<=count5; i++){
			Xo[i]=(Xn[i]*c1*c2)+(Yn[i]*(c1*s2*s3-s1*c3))+(Zn[i]*(c1*s2*c3+s1*s3));
			Yo[i]=(Xn[i]*s1*c2)+(Yn[i]*(s1*s2*s3+c1*c3))+(Zn[i]*(s1*s2*c3-c1*s3));
			Zo[i]=(Xn[i]*(-s2))+(Yn[i]*c2*s3)+(Zn[i]*c2*c3);
			dist = dist + sqrt (pow((Xo[i]-Xm[i]), 2) + pow((Yo[i]-Ym[i]), 2) + pow((Zo[i]-Zm[i]), 2) );
		}
		dist = dist/parseFloat(count5+1);
	}
}

function registration_three_nearest_neighbor_1 ()
{
	for (ii=0; ii<=cell_num_t[t+1]; ii++){
		if (cell_share1[(t+1)*c_num+ii]!=1){
			dist2 = pow((rfx[t*c_num+i]-rfx[(t+1)*c_num+ii]), 2)
			       	+ pow((rfy[t*c_num+i]-rfy[(t+1)*c_num+ii]), 2)
			       	+ pow((rfz[t*c_num+i]-rfz[(t+1)*c_num+ii]), 2);

			if (dist2<min_dist2[2]){
				min_dist2[2]=dist2;
				min_id[2] = cell_order[(t+1)*c_num+ii];

				if (dist2<min_dist2[1]){
					min_dist2[2]=min_dist2[1];
					min_dist2[1]=dist2;
					min_id[2]=min_id[1];
					min_id[1] = cell_order[(t+1)*c_num+ii];

					if (dist2<min_dist2[0]){
						min_dist2[1]=min_dist2[0];
						min_dist2[0]=dist2;
						min_id[1]=min_id[0];
						min_id[0] = cell_order[(t+1)*c_num+ii];
						min_id_order=ii;
					}
				}
			}
		}
	}			
}

function registration_three_nearest_neighbor_2 ()
{
	for (ii=0; ii<=cell_num_t[t]; ii++){
		if (cell_share1[t*c_num+ii]!=1){
			dist2 = pow((rfx[(t+1)*c_num+i]-rfx[t*c_num+ii]), 2)
			       	+ pow((rfy[(t+1)*c_num+i]-rfy[t*c_num+ii]), 2)
			       	+ pow((rfz[(t+1)*c_num+i]-rfz[t*c_num+ii]), 2);

			if (dist2<min_dist2[2]){
				min_dist2[2]=dist2;
				min_id[2] = cell_order[t*c_num+ii];

				if (dist2<min_dist2[1]){
					min_dist2[2]=min_dist2[1];
					min_dist2[1]=dist2;
					min_id[2]=min_id[1];
					min_id[1] = cell_order[t*c_num+ii];

					if (dist2<min_dist2[0]){
						min_dist2[1]=min_dist2[0];
						min_dist2[0]=dist2;
						min_id[1]=min_id[0];
						min_id[0] = cell_order[t*c_num+ii];
					}
				}
			}
		}
	}			
}

function chofuku_judge ()
{
	for (i=0; i<count7; i++){
		for (ii=i+1; ii<count7; ii++){
			if (id_chofuku[i]==id_chofuku[ii]){
				id_chofk_num[count8]=id_chofuku[i];
				count8++;
			}
		}
	}
}


function input_file_reading (input)
{
	str1 = File.openAsString(input);
	str2 = split(str1, "\t\n\r,");
	s_num = str2.length;
	p_num = s_num/5;
	print("input file name = ", input);
	print("Particle number = ", p_num);
	print("Other values in the input file: ", s_num, str2[0], str2[s_num-1]);

	count0=0;
    count1=0;
    for (i=0; i<p_num; i++)
    {
		buf_t = parseInt(str2[i*5+3]);

		if (i>0 && buf_t>end_t){
			count0++;
	       	count1=0;
		   	if (buf_t!=end_t+1){
				print("Error. Time is skipped!");
				exit(0);
			}
		}

        buf_id = parseInt(str2[i*5+4]) -1;
        end_t = buf_t;
        rfx[count0*c_num+count1] = parseFloat(str2[i*5]);
	    rfy[count0*c_num+count1] = parseFloat(str2[i*5+1]);
	    rfz[count0*c_num+count1] = parseFloat(str2[i*5+2]);;
        cell_order[count0*c_num+count1]=buf_id;
        cell_num_t[count0]=count1;
        count1++;
    }
    print("cell_number = ", cell_num_t[0], "(before)", cell_num_t[1], "(after)");
	max_cell_number=cell_num_t[0];
	if (cell_num_t[1]>max_cell_number){
		max_cell_number=cell_num_t[1];
	}
	if (max_cell_number>=c_num){	
		print("Error in cell number.");
		exit(0);
	}
}

function xyz_scaling ()
{
	//xyz coordinates are modified according to x_scaling, y_scaling, and z_scaling
	//Then, the longenst axis among x, y, and z axes is modifed to be the length = 100.
	//This scaling is important for MCMC, because the reject or accept process is affected by system scales.

	t=0;
    for (i=0; i<=cell_num_t[t]; i++)
    {
		rfx[t*c_num+i]=rfx[t*c_num+i]*x_scaling;
		rfy[t*c_num+i]=rfy[t*c_num+i]*y_scaling;
		rfz[t*c_num+i]=rfz[t*c_num+i]*z_scaling;
		if (i==0){
			max_x = rfx[t*c_num+i];
			min_x = rfx[t*c_num+i];
			max_y = rfy[t*c_num+i];
			min_y = rfy[t*c_num+i];
			max_z = rfz[t*c_num+i];
			min_z = rfz[t*c_num+i];
		}
		if (rfx[t*c_num+i]>max_x){	max_x=rfx[t*c_num+i];	}
		if (rfx[t*c_num+i]<min_x){	min_x=rfx[t*c_num+i];	}
		if (rfy[t*c_num+i]>max_y){	max_y=rfy[t*c_num+i];	}
		if (rfy[t*c_num+i]<min_y){	min_y=rfy[t*c_num+i];	}
		if (rfz[t*c_num+i]>max_z){	max_z=rfz[t*c_num+i];	}
		if (rfz[t*c_num+i]<min_z){	min_z=rfz[t*c_num+i];	}
	}
	longest_axis = (max_x-min_x);
	if ((max_y-min_y)>longest_axis){
		longest_axis = max_y-min_y;
	}
	if ((max_z-min_z)>longest_axis){
		longest_axis = max_z-min_z;
	}
	scaling_val = 100.0/longest_axis;
	print("Scaling value for xyz-coordinates = ", scaling_val);

	for (t=0; t<=1; t++)
    {
    	for (i=0; i<=cell_num_t[t]; i++)
        {
			if (t==0){
				rfx[t*c_num+i]=rfx[t*c_num+i]*scaling_val;
				rfy[t*c_num+i]=rfy[t*c_num+i]*scaling_val;
				rfz[t*c_num+i]=rfz[t*c_num+i]*scaling_val;
			}
			if (t==1){			
	        	rfx[t*c_num+i]=rfx[t*c_num+i]*x_scaling2*scaling_val;
	        	rfy[t*c_num+i]=rfy[t*c_num+i]*y_scaling2*scaling_val;
	        	rfz[t*c_num+i]=rfz[t*c_num+i]*z_scaling2*scaling_val;
			}
		}
	}

	dist_threshold = dist_thres_1*x_scaling*scaling_val;
}

function xyz_scale_reverse ()
{
	//The scales of xyz-coordinates are revesely modifed, i.e. a reversal of the xyz_scaling () function.
	//This reversal is based  not only on scaling_val but also x_scaling in the before-image. 
	
	scale_rev = 1.0/(x_scaling*scaling_val);
	for (t=0; t<=1; t++)
    {
    	for (i=0; i<=cell_num_t[t]; i++)
        {
			if (t==0){
				rfx[t*c_num+i]=rfx[t*c_num+i]*scale_rev;
				rfy[t*c_num+i]=rfy[t*c_num+i]*scale_rev;
				rfz[t*c_num+i]=rfz[t*c_num+i]*scale_rev;
			}
			if (t==1){			
	        	rfx[t*c_num+i]=rfx[t*c_num+i]*scale_rev;
	        	rfy[t*c_num+i]=rfy[t*c_num+i]*scale_rev;
	        	rfz[t*c_num+i]=rfz[t*c_num+i]*scale_rev;
			}
		}
	}
}

function landmark_particle_extraction ()
{
	for (t=0; t<=1; t++){
		for (i=0; i<=cell_num_t[t]; i++){
		       	cell_share1[t*c_num+i]=0;
		}
	}
	t=0;
	count2=-1;
	for (i=0; i<=cell_num_t[t]; i++){
		buf_id=cell_order[t*c_num+i];
		for (j=0; j<=cell_num_t[t+1]; j++){
			if (buf_id==cell_order[(t+1)*c_num+j]){
				cell_share1[t*c_num+i]=1;
				cell_share1[(t+1)*c_num+j]=1;
				cell_pair[i]=j;
				count2++;
			}
		}
	}
	landmark_num = count2 + 1;
	print("Number of landmark particles = ", landmark_num);
}

function translational_motion_removal ()
{
	t=1;

	//centroid of landmark particles in before-image
	count0=0;
	centroid_x1=0.0;
    centroid_y1=0.0;
    centroid_z1=0.0;
	for (i=0; i<=cell_num_t[t-1]; i++){
		if (cell_share1[(t-1)*c_num+i]==1){	//choose landmark particles
			count0++;
			centroid_x1=centroid_x1+rfx[(t-1)*c_num+i];
			centroid_y1=centroid_y1+rfy[(t-1)*c_num+i];
			centroid_z1=centroid_z1+rfz[(t-1)*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x1=centroid_x1/parseFloat(count0);
		centroid_y1=centroid_y1/parseFloat(count0);
		centroid_z1=centroid_z1/parseFloat(count0);
	}

	//centroid of landmark particles in after-image
	count0=0;
	centroid_x2=0.0;
    centroid_y2=0.0;
    centroid_z2=0.0;
	for (i=0; i<=cell_num_t[t]; i++){
		if (cell_share1[t*c_num+i]==1){	//choose landmark particles
			count0++;
			centroid_x2=centroid_x2+rfx[t*c_num+i];
			centroid_y2=centroid_y2+rfy[t*c_num+i];
			centroid_z2=centroid_z2+rfz[t*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x2=centroid_x2/parseFloat(count0);
		centroid_y2=centroid_y2/parseFloat(count0);
		centroid_z2=centroid_z2/parseFloat(count0);
	}

	//modification of xyz coordinates in after-image according to the centroid translation between before-image and after-image
	cent_displace_x=centroid_x2-centroid_x1;
	cent_displace_y=centroid_y2-centroid_y1;
	cent_displace_z=centroid_z2-centroid_z1;
	for (i=0; i<=cell_num_t[t]; i++){
		rfx[t*c_num+i]=rfx[t*c_num+i]-cent_displace_x;
		rfy[t*c_num+i]=rfy[t*c_num+i]-cent_displace_y;
		rfz[t*c_num+i]=rfz[t*c_num+i]-cent_displace_z;
	}
}

function translation_movement_to_origin ()
{
	t=1;

	//centroid of landmark particles in before-image
	count0=0;
    centroid_x1=0.0;
    centroid_y1=0.0;
    centroid_z1=0.0;
	for (i=0; i<=cell_num_t[t-1]; i++){	//choose landmark particles
		if (cell_share1[(t-1)*c_num+i]==1){
			count0++;
			centroid_x1=centroid_x1+rfx[(t-1)*c_num+i];
			centroid_y1=centroid_y1+rfy[(t-1)*c_num+i];
			centroid_z1=centroid_z1+rfz[(t-1)*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x1=centroid_x1/parseFloat(count0);
		centroid_y1=centroid_y1/parseFloat(count0);
		centroid_z1=centroid_z1/parseFloat(count0);
		count4=-1;
		//xyz coordinates are moved to (0,0,0).
		for (i=0; i<=cell_num_t[t-1]; i++){
			if (cell_share1[(t-1)*c_num+i]==1){
				count4++;
				Xm[count4]=rfx[(t-1)*c_num+i]-centroid_x1;
			    Ym[count4]=rfy[(t-1)*c_num+i]-centroid_y1;
			    Zm[count4]=rfz[(t-1)*c_num+i]-centroid_z1;
			}				
		}
	}

	//centroid of landmark particles in after-image
	count0=0;
	centroid_x2=0.0;
    centroid_y2=0.0;
    centroid_z2=0.0;
	for (i=0; i<=cell_num_t[t]; i++){
		if (cell_share1[t*c_num+i]==1){	//choose landmark particles
			count0++;
			centroid_x2=centroid_x2+rfx[t*c_num+i];
			centroid_y2=centroid_y2+rfy[t*c_num+i];
			centroid_z2=centroid_z2+rfz[t*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x2=centroid_x2/parseFloat(count0);
		centroid_y2=centroid_y2/parseFloat(count0);
		centroid_z2=centroid_z2/parseFloat(count0);
		count4=-1;
		//xyz coordinates are moved to (0,0,0).
		for (i=0; i<=cell_num_t[t]; i++){
			if (cell_share1[t*c_num+i]==1){
				count4++;
				Xn[count4]=rfx[t*c_num+i]-centroid_x2;
			    Yn[count4]=rfy[t*c_num+i]-centroid_y2;
			   	Zn[count4]=rfz[t*c_num+i]-centroid_z2;
			}				
		}
		count5=count4;		//number of landmark particles
	}
}

function output_conditions (input_name)
{
	output_7 = dir_input+"output_condition_for_registration.txt";

	File.saveString("input_file_name = "+input_name+"\n", output_7);

	File.append("xyz scaling for before-image = "+toString(x_scaling)+", "+toString(y_scaling)+", "+toString(z_scaling), output_7);
	File.append("xyz scaling for after-image = "+toString(x_scaling2)+", "+toString(y_scaling2)+", "+toString(z_scaling2), output_7);
	File.append("xyz scaling for all image = "+toString(scaling_val), output_7);
	File.append("xyz longest axis length after scaling = "+toString(longest_axis*scaling_val)+"\n", output_7); 
	File.append("number of landmark particles = "+toString(landmark_num)+"\n", output_7);
	File.append("****** MCMC parameters ******", output_7);
	File.append("	iteration = "+toString(end_iteration)+"\n"+"	kBT="+toString(kBT)+"\n"+"	mag_xyz="+toString(mag_xyz)+"\n"+"	uniform random number="+toString(mag)+"\n", output_7);
	File.append("	segment number of angles = "+toString(d_ang), output_7);
	File.append("		which means "+toString(d_ang*d_ang*d_ang)+" initial conditions tried\n", output_7);
	File.append("	threshold value as mean distance of paired particles for MCMC", output_7);
	File.append("		= "+toString(dist_thres_1)+"(original scale)\n"+"		= "+toString(dist_threshold)+"(revised scale)\n", output_7);

	//*** save log file ***
	selectWindow("Log");
	saveAs("Text", dir_input+"/Log_registration.txt");		
}

function output_xyz_coordinate_after_rotation ()
{
	output_5 = dir_input+"output_registered_particles_landmarks.txt";
	output_6 = dir_input+"output_registered_particles_all.txt";

	File.saveString("before-image, after-image, distances, xyz-coordinates (before & after)\n", output_5);
	File.append("id0, id1, dist, x0, y0, z0, x1, y1, z1", output_5);
	File.saveString("before-image, after-image, distances, xyz-coordinates (before & after)\n", output_6);
	File.append("id0, id1, dist, x0, y0, z0, x1, y1, z1", output_6);

	t=0;
	for (i=0; i<=cell_num_t[t]; i++){
		ii=cell_pair[i];
		if (cell_share1[t*c_num+i]==1){
			buf_id=cell_order[(t+1)*c_num+ii];
		}
		else{
			buf_id=paired_id[i];
		}
		dist2 = pow((rfx[t*c_num+i]-rfx[(t+1)*c_num+ii]), 2)
		       	+ pow((rfy[t*c_num+i]-rfy[(t+1)*c_num+ii]), 2)
		       	+ pow((rfz[t*c_num+i]-rfz[(t+1)*c_num+ii]), 2);
		//all particles
		str3 = toString(cell_order[t*c_num+i]+1)+", "+toString(cell_order[(t+1)*c_num+ii]+1)+", "+toString(sqrt(dist2))+", "+toString(rfx[t*c_num+i])+", "+toString(rfy[t*c_num+i])+", "+toString(rfz[t*c_num+i])+", "+toString(rfx[(t+1)*c_num+ii])+", "+toString(rfy[(t+1)*c_num+ii])+", "+toString(rfz[(t+1)*c_num+ii]);
		File.append(str3, output_6);

		//landmark particles only
		if (cell_share1[t*c_num+i]==1){
			str3 = toString(cell_order[t*c_num+i]+1)+", "+toString(cell_order[(t+1)*c_num+ii]+1)+", "+toString(sqrt(dist2))+", "+toString(rfx[t*c_num+i])+", "+toString(rfy[t*c_num+i])+", "+toString(rfz[t*c_num+i])+", "+toString(rfx[(t+1)*c_num+ii])+", "+toString(rfy[(t+1)*c_num+ii])+", "+toString(rfz[(t+1)*c_num+ii]);
			File.append(str3, output_5);
		}
	}
}

function output_pairing_particles ()
{
	output_3 = dir_input+"output_registered_particles_3neighbors.txt";
	output_4 = dir_input+"output_registered_particles_1neighbors.txt";

	t=0;

	//before-image (i.e. t=0) is the reference.
	File.saveString("before-image, after-image (3-neighbors), distances\n", output_3);
	File.append("id0, id1, id2, id3, dist1, dist2, dist3", output_3);
	File.saveString("before-image, after-image (1-neighbor), distances\n", output_4);
	File.append("id0, id1, dist1", output_4);
	count7=0;
	for (i=0; i<=cell_num_t[t]; i++){
		min_dist2[0]=10000000;
		min_dist2[1]=100000000;
		min_dist2[2]=1000000000;
		min_id[0]=10000;
		min_id[1]=10000;
		min_id[2]=10000;
		if (cell_share1[t*c_num+i]!=1){
			registration_three_nearest_neighbor_1 ();
			str3 = toString(cell_order[t*c_num+i]+1)+", "+toString(min_id[0]+1)+", "+toString(min_id[1]+1)+", "+toString(min_id[2]+1)+", "+toString(sqrt(min_dist2[0]))+", "+toString(sqrt(min_dist2[1]))+", "+toString(sqrt(min_dist2[2]));
			File.append(str3, output_3);
			str3 = toString(cell_order[t*c_num+i]+1)+", "+toString(min_id[0]+1)+", "+toString(sqrt(min_dist2[0]));
			File.append(str3, output_4);

			id_chofuku[count7]=min_id[0];
			paired_id[i]=min_id[0];
			cell_pair[i]=min_id_order;
			count7++;
		}
	}
	
	//Searching for multiplly counted particles
	count8=0;
	chofuku_judge ();
	str3 = "multiplly counted particles = "+toString(count8);
	File.append(str3, output_3);
	File.append(str3, output_4);
	for (i=0; i<count8; i++){
		if (i==0){
			str3 = toString(id_chofk_num[i]+1);
		}
		else{
			str3 = str3+" "+toString(id_chofk_num[i]+1);
		}
	}
	File.append(str3, output_3);
	File.append(str3, output_4);

	//after-image (i.e. t=1) is the reference.
	File.append("\nafter-image, before-image (3-neighbors), distances", output_3);
	File.append("id0, id1, id2, id3, dist1, dist2, dist3", output_3);
	File.append("\nafter-image, before-image (1-neighbor), distances", output_4);
	File.append("id0, id1, dist1", output_4);
	count7=0;
	for (i=0; i<=cell_num_t[t+1]; i++){
		min_dist2[0]=10000000;
		min_dist2[1]=100000000;
		min_dist2[2]=1000000000;
		min_id[0]=10000;
		min_id[1]=10000;
		min_id[2]=10000;
		if (cell_share1[(t+1)*c_num+i]!=1){
			registration_three_nearest_neighbor_2 ();
			str3 = toString(cell_order[(t+1)*c_num+i]+1)+", "+toString(min_id[0]+1)+", "+toString(min_id[1]+1)+", "+toString(min_id[2]+1)+", "+toString(sqrt(min_dist2[0]))+", "+toString(sqrt(min_dist2[1]))+", "+toString(sqrt(min_dist2[2]));
			File.append(str3, output_3);
			str3 = toString(cell_order[(t+1)*c_num+i]+1)+", "+toString(min_id[0]+1)+", "+toString(sqrt(min_dist2[0]));
			File.append(str3, output_4);
			id_chofuku[count7]=min_id[0];
			count7++;
		}
	}

	//Searching for multiplly counted particles
	count8=0;
	chofuku_judge ();
	str3 = "multiplly counted particles = "+toString(count8);
	File.append(str3, output_3);
	File.append(str3, output_4);
	for (i=0; i<count8; i++){
		if (i==0){
			str3 = toString(id_chofk_num[i]+1);
		}
		else{
			str3 = str3+" "+toString(id_chofk_num[i]+1);
		}
	}
	File.append(str3, output_3);
	File.append(str3, output_4);
}




