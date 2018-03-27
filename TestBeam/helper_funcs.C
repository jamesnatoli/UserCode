#include "slimanalysis.h"
#include <string>

/* Documentation */
// Fills the array thetas with the appropriate angle to use for each tile
void fill_Rot_Array() {
  for (int chan = 0; chan < NUMCHAN; chan++) {
    thetas[chan] = calc_theta(chan); // fill the theta array
    for (int pt = 0; pt < 4; pt++) {
      // Fill the Rotated Arrays
      rot_fiducialX[chan][pt] = rotate_Point(fiducialX[chan][pt], fiducialY[chan][pt]
					     , chan, 'x');
      rot_fiducialY[chan][pt] = rotate_Point(fiducialX[chan][pt], fiducialY[chan][pt]
					     , chan, 'y');
    } // points
  } // channels 
}

double calc_theta(int channel_num) {
  // 200 => BL                                                  
  if ( channel_num == 2)
    return atan(abs(fiducialY[channel_num][1] - fiducialY[channel_num][0]) /
                abs(fiducialX[channel_num][1] - fiducialX[channel_num][0]));
  // 260, Fingers1, 2, 3, 81S => TR
  else if (channel_num == 0 || channel_num == 3 || channel_num == 4
           || channel_num == 5 || channel_num == 7)
    return atan(abs(fiducialY[channel_num][3] - fiducialY[channel_num][2]) /
                abs(fiducialX[channel_num][3] - fiducialX[channel_num][2]));
  // 260-2p and Finger 4 => TL
  else if (channel_num == 1 || channel_num == 6)
    return atan(abs(fiducialY[channel_num][3] - fiducialY[channel_num][2]) /
                abs(fiducialX[channel_num][3] - fiducialX[channel_num][2]));
  return 0;
}

bool isFiducial(int i, float x_hit, float y_hit) {
  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((fiducialY[i][k]< y_hit && fiducialY[i][j]>=y_hit) ||
         (fiducialY[i][j]< y_hit && fiducialY[i][k]>=y_hit)) &&
        (fiducialX[i][k]<=x_hit || fiducialX[i][j]<=x_hit)) {
      oddNodes^=(fiducialX[i][k]+(y_hit-fiducialY[i][k])/
                 (fiducialY[i][j]-fiducialY[i][k])*
                 (fiducialX[i][j]-fiducialX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}

bool isRotFiducial(int i, float x_hit, float y_hit) {
  if (rot_fiducialX[0][0] == 0 || rot_fiducialY[0][0] == 0)
    cout << "DANGER WILL ROBINSON" << endl;

  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((rot_fiducialY[i][k]< y_hit && rot_fiducialY[i][j]>=y_hit) ||
         (rot_fiducialY[i][j]< y_hit && rot_fiducialY[i][k]>=y_hit)) &&
        (rot_fiducialX[i][k]<=x_hit || rot_fiducialX[i][j]<=x_hit)) {
      oddNodes^=(rot_fiducialX[i][k]+(y_hit-rot_fiducialY[i][k])/
                 (rot_fiducialY[i][j]-rot_fiducialY[i][k])*
                 (rot_fiducialX[i][j]-rot_fiducialX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}
/* 
bool isFiducial(int i, float x_hit, float y_hit, float arrX[NUMCHAN][4], float arrY[NUMCHAN][4]) {
  // Find if it is fiducial!
  int polyCorners = 4;
  int k, j=polyCorners-1 ;
  bool oddNodes=false;
  for (k=0; k<polyCorners; k++) {
    if (((arrY[i][k]< y_hit && arrY[i][j]>=y_hit) ||
         (arrY[i][j]< y_hit && arrY[i][k]>=y_hit)) &&
        (arrX[i][k]<=x_hit || arrX[i][j]<=x_hit)) {
      oddNodes^=(arrX[i][k]+(y_hit-arrY[i][k])/
                 (arrY[i][j]-arrY[i][k])*
                 (arrX[i][j]-arrX[i][k])<x_hit);
    }
    j=k;
  }
  return oddNodes;
}
*/

/* double rotate_Point(double point_X, double point_Y, double theta, char xy) { */
/*   if (toupper(xy) == 'X') { */
/*     return ((point_X*cos(theta)) + (point_Y*sin(theta))); */
/*   } */
/*   else { */
/*     return ((-point_X*sin(theta)) + (point_Y*cos(theta))); */
/*   } */
/* } */

double rotate_Point(double point_X, double point_Y, int channel_num, char xy) {
  if (toupper(xy) == 'X') {
    return ((point_X*cos(thetas[channel_num])) + (point_Y*sin(thetas[channel_num])));
  }
  else {
    return ((-point_X*sin(thetas[channel_num])) + (point_Y*cos(thetas[channel_num])));
  }
}

/*
bool is_Left_Or_Right( int chan, double px, double py) {
  if (px < min(rot_fiducialX[chan][0], rot_fiducialX[chan][2])
      || px > max(rot_fiducialX[chan][1], rot_fiducialX[chan][3]))
    return false;
  else if (px > max(rot_fiducialX[chan][0], rot_fiducialX[chan][2])
	   || px > min(rot_fiducialX[chan][1], rot_fiducialX[chan][3]))
    return true;
  else {
    double m = (abs(rot_fiducialY[chan][2] - rot_fiducialY[chan][0])
		/ abs(rot_fiducialX[chan][2] - rot_fiducialX[chan][0]));
    if ((m * (px - rot_fiducialX[chan][0])

	 } */

