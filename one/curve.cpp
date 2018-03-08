#include "curve.h"
#include "extra.h"
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}
    

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // TODO:
    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Make sure that this function computes all the appropriate
    // Vector3fs for each CurvePoint: V,T,N,B.
    // [NBT] should be unit and orthogonal.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

	// preallocate curve with (|P|/4)*steps points
	Curve C((P.size()/4) * steps);

	//TODO: do this for each 4 points successivly where last point is first point of next 
	for (unsigned i = 0; i < P.size(); i += 4) {
		//TODO: bezier EQ ( 4 points )
		// p = (1-t)^3 *P0 + 3*t*(1-t)^2*P1 + 3*t^2*(1-t)*P2 + t^3*P3 

		double stepsize = 1 / (double) steps;
		double t = 0;

		for (unsigned j = 0; j < steps; j++) {

            CurvePoint cp;
			// cubic bezier equation
			Vector3f pp = pow((1 - t), 3)*P[i] + 3 * t*pow((1 - t), 2)*P[i + 1] + 3 *t*t*(1 - t)*P[i + 2] + pow(t, 3)*P[i + 3];
			cp.V = pp;

			// Tangent vector is first derivative
			cp.T = Vector3f(0, 0, 0);
			// Normal vector is second derivative
			cp.N = Vector3f(0, 0, 0);
			// Finally, binormal is facing up.
			cp.B = Vector3f(0, 0, 1);

            // add point to end of curve
            C.push_back(cp);
            
			t += stepsize;
		}
	}

    return C;

}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    // TODO:
    // It is suggested that you implement this function by changing
    // basis from B-spline to Bezier.  That way, you can just call
    // your evalBezier function.

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

    Matrix4f bezierBasis (
        1, -3,  3, -1,
        0,  3, -6,  3,
        0,  0,  3, -3,
        0,  0,  0,  1
    );

    Matrix4f bsplineBasis (
        1/6.0f, -3/6.0f,  3/6.0f, -1/6.0f,
        4/6.0f,       0, -6/6.0f,  3/6.0f,
        1/6.0f,  3/6.0f,  3/6.0f, -3/6.0f, 
             0,       0,       0,  1/6.0f
    );

    // conversion matrix is B1 * B2^-1
    Matrix4f conversionMatrix = bsplineBasis * bezierBasis.inverse();
    vector<Vector3f> newControlPoints;

    for(unsigned i=0; i< P.size() - 3; i++) {
        // for each row (create new control point)
        for(unsigned j=0; j < 4; j++) {
            Vector3f newControlPoint(0.0f);
            // for each col
            for(unsigned k=0; k < 4; k++) {
                newControlPoint += conversionMatrix[(j*4) + k] * P[i+k];
            }
            newControlPoints.push_back(newControlPoint);
        }
    }

    newControlPoints.pop_back();

    cerr <<  "NUMBER OF CONTROL POINTS INTO BEZIER" << newControlPoints.size() << "\n\n";

    return evalBezier(newControlPoints, steps);
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i ) {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize ) {
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

