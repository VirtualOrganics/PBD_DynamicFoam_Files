#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "GL/glut.h"
#include "Demos/Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "GenericConstraintsModel.h"
#include <iostream>
#include "Demos/Simulation/TimeStepController.h"
#include "Demos/Visualization/Visualization.h"
#include "Demos/Utils/Utilities.h"
#include "Demos/Simulation/DistanceFieldCollisionDetection.h"
#include "Demos/Utils/OBJLoader.h"
#include "Demos/Utils/Timing.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
	#define new DEBUG_NEW 
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;

void timeStep ();
void buildModel ();
void createMesh();
void render ();
void cleanup();
void reset();
void initShader();
void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end);
void TW_CALL setTimeStep(const void *value, void *clientData);
void TW_CALL getTimeStep(void *value, void *clientData);
void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData);
void TW_CALL getVelocityUpdateMethod(void *value, void *clientData);
void TW_CALL setStiffness(const void *value, void *clientData);
void TW_CALL getStiffness(void *value, void *clientData);
void TW_CALL setBendingStiffness(const void *value, void *clientData);
void TW_CALL getBendingStiffness(void *value, void *clientData);

GenericConstraintsModel model;
TimeStepController simulation;
DistanceFieldCollisionDetection dfcd;

/*
const int nRows = 30;
const int nCols = 30;
*/
const int nRows = 3;
const int nCols = 3;
const Real width = 4.0;
const Real height = 4.0;
bool doPause = true;
std::vector<unsigned int> selectedParticles;
std::vector<unsigned int> currentVertices;
Vector3r oldMousePos;
Shader *shader;
string exePath;
string dataPath;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS

	exePath = Utilities::getFilePath(argv[0]);
	// dataPath = exePath + "/" + std::string(PBD_DATA_PATH);
	dataPath = exePath + std::string(PBD_DATA_PATH);

	// OpenGL
	MiniGL::init (argc, argv, 1024, 768, 0, 0, "Michel Van de gaer");
	MiniGL::initLights ();
	MiniGL::initTexture();
	MiniGL::setClientIdleFunc (50, timeStep);		
	MiniGL::setKeyFunc(0, 'r', reset);
	MiniGL::setSelectionFunc(selection);
	initShader();

	buildModel ();

	MiniGL::setClientSceneFunc(render);			
	//MiniGL::setViewport (40.0f, 0.1f, 500.0f, Vector3r (5.0, 10.0, 30.0), Vector3r (5.0, 0.0, 0.0));

	MiniGL::setViewport (40.0f, 0.1f, 500.0f, Vector3r (0.0, 40.0, 0.0), Vector3r (0.0, 0.0, 0.0));

	// MiniGL::setViewport (40.0f, 0.1f, 500.0f, Vector3r (0.0, 0.0, 0.0), Vector3r (0.0, 0.0, 0.0));

	TwAddVarRW(MiniGL::getTweakBar(), "Pause", TW_TYPE_BOOLCPP, &doPause, " label='Pause' group=Simulation key=SPACE ");
	TwAddVarCB(MiniGL::getTweakBar(), "TimeStepSize", TW_TYPE_REAL, setTimeStep, getTimeStep, &model, " label='Time step size'  min=0.0 max = 1.0 step=0.001 precision=6 group=Simulation ");
	TwType enumType = TwDefineEnum("VelocityUpdateMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "VelocityUpdateMethod", enumType, setVelocityUpdateMethod, getVelocityUpdateMethod, &simulation, " label='Velocity update method' enum='0 {First Order Update}, 1 {Second Order Update}' group=Simulation");
	TwAddVarCB(MiniGL::getTweakBar(), "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, &model, " label='Stiffness'  min=0.0 step=0.1 precision=6 group='Distance constraints' ");
	TwAddVarCB(MiniGL::getTweakBar(), "BendingStiffness", TW_TYPE_REAL, setBendingStiffness, getBendingStiffness, &model, " label='Bending stiffness'  min=0.0 step=0.01 precision=6 group=Bending ");

	glutMainLoop ();	

	cleanup ();
	
	Timing::printAverageTimes();

	return 0;
}

void initShader()
{
	std::string vertFile = dataPath + "/shaders/vs_smoothTex.glsl";
	std::string fragFile = dataPath + "/shaders/fs_smoothTex.glsl";
	shader = MiniGL::createShader(vertFile, "", fragFile);

	if (shader == NULL) {
		// std::cout << "shader is NULL" << "\n";
		return;
	}

	shader->begin();
	shader->addUniform("modelview_matrix");
	shader->addUniform("projection_matrix");
	shader->addUniform("surface_color");
	shader->addUniform("shininess");
	shader->addUniform("specular_factor");
	shader->end();
}

void cleanup()
{
	delete TimeManager::getCurrent();
	delete shader;
}

void reset()
{
	Timing::printAverageTimes();
	Timing::reset();

	model.reset();
	simulation.reset();
	TimeManager::getCurrent()->setTime(0.0);
}

void mouseMove(int x, int y)
{
	Vector3r mousePos;
	MiniGL::unproject(x, y, mousePos);
	const Vector3r diff = mousePos - oldMousePos;

	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	ParticleData &pd = model.getParticles();
	for (unsigned int j = 0; j < selectedParticles.size(); j++)
	{
		pd.getVelocity(selectedParticles[j]) += 5.0*diff/h;
	}
	oldMousePos = mousePos;
	// fprintf (stderr,"%d %d\n",x,y);
}

void selection(const Eigen::Vector2i &start, const Eigen::Vector2i &end)
{
	Eigen::Vector2i ss, ee;
	std::vector<unsigned int> hits;
	selectedParticles.clear();
	currentVertices.clear();
	ParticleData &pd = model.getParticles();

	ss[0] = ss[1] = 0;
	ee[0] = 1024;
	ee[1] = 768;
	// Selection::selectRect(start, end, &pd.getPosition(0), &pd.getPosition(pd.size() - 1), selectedParticles);
	Selection::selectRect(ss, ee, &pd.getPosition(0), &pd.getPosition(pd.size() - 1), selectedParticles, currentVertices);
	if (selectedParticles.size() > 0)
		MiniGL::setMouseMoveFunc(GLUT_MIDDLE_BUTTON, mouseMove);
	else
		MiniGL::setMouseMoveFunc(-1, NULL);

	// MiniGL::unproject(end[0], end[1], oldMousePos);
	MiniGL::unproject(ee[0], ee[1], oldMousePos);

	// fprintf (stderr,"selection: %d %d\n",selectedParticles.size(),currentVertices.size());

	// fprintf (stderr,"%d %d %d %d\n",start[0],start[1],end[0],end[1]);
	// fprintf (stderr,"%d %d %d %d\n",ss[0],ss[1],ee[0],ee[1]);
}

Vector3r calc_int (Vector3r obj1, Vector3r obj2, Real rad1, Real rad2, Real r0, int which)
{
	Vector3r	aa, bb, cc, dd;	
	Real		ee, ff;

	// calculate the intersection(s) of two circles
	aa[0] = 0.5 * (obj1[0]+obj2[0]);
	aa[2] = 0.5 * (obj1[2]+obj2[2]);

	bb[0] = ((rad1*rad1) - (rad2*rad2)) / (2*r0*r0) * (obj2[0]-obj1[0]);
	bb[2] = ((rad1*rad1) - (rad2*rad2)) / (2*r0*r0) * (obj2[2]-obj1[2]);

	ff = (rad1*rad1)-(rad2*rad2);
	ee = 0.5 * sqrt ((2 * ((rad1*rad1)+(rad2*rad2))/(r0*r0)) -
			((ff*ff)/(r0*r0*r0*r0)) - 1);

	cc[0] = ee * (obj2[2] - obj1[2]);
	cc[2] = ee * (obj1[0] - obj2[0]);

	if (which == 0) {
		dd[0] = aa[0] + bb[0] + cc[0];
		dd[2] = aa[2] + bb[2] + cc[2];
	} else {
		dd[0] = aa[0] + bb[0] - cc[0];
		dd[2] = aa[2] + bb[2] - cc[2];
	}

	return (dd);
}

void timeStep ()
{
	if (doPause)
		return;

	ParticleData &pd = model.getParticles();
	CurrentData &cd = model.getCurrents();
	Real	r0, r1, rad0, rad1, rad2;
	Real	nInt;
	Real	tIme;
	Real	aa, bb;
	Real	xx, yy;
	Vector3r    pos0, pos1, pos2, pos_a, pos_b, pos_a1, pos_b1;
	Vector3r	pos_x1, pos_x2;
	Vector3r    pos, pos3, diff;
	unsigned int	i, k, m, n, a, b, cc;
	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	float grn[4] = { 0.0f, 0.8f, 0.0f, 1 };
	float blu[4] = { 0.0f, 0.0f, 0.8f, 1 };
	float blk[4] = { 0.0f, 0.0f, 0.0f, 1 };

	// Simulation code
	for (i = 0; i < 8; i++) {
	// for (unsigned int i = 0; i < 1; i++) {
		// fprintf (stderr,"%d\n",i);
		simulation.step(model);
	}

	// calculate the intersections
	for (unsigned int j = 0; j < selectedParticles.size(); j++) {

		pd.setnInt (selectedParticles[j], 0);

		pos0 = pd.getPosition(selectedParticles[j]);
		rad0 = pd.getRadius(selectedParticles[j]);

		for (i = 0 ; i < j ; i++) {

			pos1 = pd.getPosition(selectedParticles[i]);
			diff = pos1 - pos0;
			r0 = sqrt ((diff[0]*diff[0])+(diff[2]*diff[2]));
	
			rad1 = pd.getRadius(selectedParticles[i]);
	
			if (r0 <= (float)(rad0+rad1)) {
				// fprintf (stderr,"%f %d %d %f %f\n",TimeManager::getCurrent ()->getTime (),j,i,r0,(float)(rad0+rad1));
				pos = calc_int (pos0, pos1, rad0, rad1, r0, 0);
				pos3 = calc_int (pos0, pos1, rad0, rad1, r0, 1);

				nInt = pd.getnInt (selectedParticles[j]);
				// fprintf (stderr,"%f\n",nInt);
				pd.setInt (selectedParticles[j], pos, (unsigned int) nInt);
				pd.addInt (selectedParticles[j], i, (unsigned int) (nInt / 2));
				nInt = pd.getnInt (selectedParticles[j]);
				pd.setInt (selectedParticles[j], pos3, (unsigned int) nInt);
				pd.addInt (selectedParticles[j], i, (unsigned int) (nInt / 2));

				nInt = pd.getnInt (selectedParticles[i]);
				pd.setInt (selectedParticles[i], pos3, (unsigned int) nInt);
				pd.addInt (selectedParticles[i], j, (unsigned int) (nInt / 2));
				nInt = pd.getnInt (selectedParticles[i]);
				pd.setInt (selectedParticles[i], pos, (unsigned int) nInt);
				pd.addInt (selectedParticles[i], j, (unsigned int) (nInt / 2));
			}
		}
	}

	/*
	for (unsigned int j = 0; j < selectedParticles.size(); j++) {
		pos0 = pd.getPosition(selectedParticles[j]);
		fprintf (stderr,"%d %f %f\n",j,pos0[0],pos0[2]);
	}
	*/


	// calculate the current network
	cd.release();
	m = 0;
	n = 0;
	for (unsigned int j = 0; j < selectedParticles.size(); j++) {
		nInt = pd.getnInt (selectedParticles[j]);

		fprintf (stderr, "%d %d - ",j,(unsigned int) nInt/2);
		a = (unsigned int) (nInt/2);
		for (k = 0 ; k < a ; k++) {
			rad0 = pd.getWhich (j, k % a);
			rad1 = pd.getWhich (j, (k+1) % a);
		/*
		for (k = 0 ; k < a-1 ; k++) {
			rad0 = pd.getWhich (j, k);
			rad1 = pd.getWhich (j, k+1);
		*/
			fprintf (stderr,"[%d %d] ",(unsigned int) rad0, (unsigned int) rad1);

			b = (unsigned int) pd.getnInt ((unsigned int) rad0);
			if (b > 0) {
				xx = yy = 0.0;
				for (cc = 0 ; cc < b/2 ; cc++) {
					rad2 = pd.getWhich ((unsigned int) rad0, cc);
					if (rad2 == rad1) {
						fprintf (stderr,"%d %d %d ",j,(unsigned int) rad0,(unsigned int) rad1);
						pos1 = pd.getPosition (selectedParticles[j]);
						xx += pos1[0];
						yy += pos1[2];
						pos1 = pd.getPosition (selectedParticles[(unsigned int) rad0]);
						xx += pos1[0];
						yy += pos1[2];
						pos1 = pd.getPosition (selectedParticles[(unsigned int) rad1]);
						xx += pos1[0];
						yy += pos1[2];
						fprintf (stderr,"yes ");
						pos0[0] = xx / 3.0;
						pos0[1] = pos0[1];
						pos0[2] = yy / 3.0;
						cd.setCurrA (n, pos0);
						n++;
					}
				}
			}
		}
		fprintf (stderr,"\n%d\n",n);

		/*
			// now determine which one of the circle intersections (points) are furthest away from the centroid just calculated.
			// These will determine the other end of a vector connecting the two points.
			for (k = 0 ; k < (unsigned int ) nInt ; k+=2) {
				// fprintf (stderr,"+ %f %f %d %d %f %f %f\n",tIme,nInt,(unsigned int)(nInt/2),k/2,rad0,pos1[0],pos1[2]);

				// get one intersection
				pos_a = pd.getInt (j, k);

				// get the other intersection
				pos_b = pd.getInt (j, k+1);

				aa = (pos_a[0]-pos[0]);
				bb = (pos_a[2]-pos[2]);
				r0 = sqrt ((aa*aa)+(bb*bb));
				aa = (pos_b[0]-pos[0]);
				bb = (pos_b[2]-pos[2]);
				r1 = sqrt ((aa*aa)+(bb*bb));
				if (r0 > r1) {
					pos2 = pos_a;
				} else {
					pos2 = pos_b;
				}
				cd.setCurrA (n, pos);
				cd.setCurrB (n, pos2);
				n++;
			}
			m++;
		}
		*/
	}
	/*
	cd.setCurrA (n, pos_x1);
	cd.setCurrB (n, pos_x2);
	n++;
	*/

	// modify particles zero's radius randomly
	/*
	rad0 = pd.getRadius(selectedParticles[0]);
	// fprintf (stderr,"getRadius: %f\n",rad0);
	rad1 = (Real) rand() / RAND_MAX;
	if (rad1 <= (Real) 0.5) {
		rad0 += (Real) 0.05;
	} else {
		rad0 -= (Real) 0.05;
	}
	if (rad0 < 1.0) {
		rad0 = 1.0;
	}
	if (rad0 > 4.0) {
		rad0 = 4.0;
	}
	pd.setRadius(selectedParticles[0], rad0);
	*/

	/*
	if (rad0 < 1.0) {
		rad0 = 4.0;
	}
	if (rad0 > 4.0) {
		rad0 = 1.0;
	}
	*/

	// modify all particle radii randomly
	for (unsigned int j = 0; j < selectedParticles.size(); j++) {
		rad0 = pd.getRadius(selectedParticles[j]);
		rad1 = (Real) rand() / RAND_MAX;
		if (rad1 <= (Real) 0.5) {
			rad0 += (Real) 0.05;
		} else {
			rad0 -= (Real) 0.05;
		}
		if (rad0 < 1.0) {
			rad0 = 1.0;
		}
		if (rad0 > 4.0) {
			rad0 = 4.0;
		}
		pd.setRadius(selectedParticles[j], rad0);
		// pd.setnInt (selectedParticles[j], 0);
	}
}

void buildModel ()
{
	TimeManager::getCurrent ()->setTimeStepSize (0.005);

	createMesh();

	simulation.setCollisionDetection(model, &dfcd);
	// cd.setTolerance(0.05);
	dfcd.setTolerance(0.00);
}

void renderTriangleModels()
{
	// Draw simulation model

	const ParticleData &pd = model.getParticles();
	float surfaceColor[4] = { 0.2f, 0.5f, 1.0f, 1 };

	if (shader)
	{
		shader->begin();
		glUniform3fv(shader->getUniform("surface_color"), 1, surfaceColor);
		glUniform1f(shader->getUniform("shininess"), 5.0f);
		glUniform1f(shader->getUniform("specular_factor"), 0.2f);

		GLfloat matrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
		glUniformMatrix4fv(shader->getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
		GLfloat pmatrix[16];
		glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
		glUniformMatrix4fv(shader->getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);
	}

	for (unsigned int i = 0; i < model.getTriangleModels().size(); i++)
	{
		// mesh 
		TriangleModel *triModel = model.getTriangleModels()[i];
		const IndexedFaceMesh &mesh = triModel->getParticleMesh();
		// Visualization::drawTexturedMesh(pd, mesh, triModel->getIndexOffset(), surfaceColor);
	}
	if (shader)
		shader->end();
}

void render ()
{
	Vector3r	pos,pos3;
	unsigned int	i, j, k;
	MiniGL::coordinateSystem();
	
	/*
	if (doPause)
		return;
	*/

	renderTriangleModels();

	float red[4] = { 0.8f, 0.0f, 0.0f, 1 };
	float grn[4] = { 0.0f, 0.8f, 0.0f, 1 };
	float blu[4] = { 0.0f, 0.0f, 0.8f, 1 };
	float blk[4] = { 0.0f, 0.0f, 0.0f, 1 };
	// const ParticleData &pd = model.getParticles();
	ParticleData &pd = model.getParticles();
	CurrentData &cd = model.getCurrents();
	Vector3r    pos0, pos1, pos2, pos_a, pos_b, pos_a1, pos_b1;
	Vector3r	pos_x1, pos_x2;
	Vector3r	diff;
	float	r0, r1, m;
	Real	rad0, rad1, rad00;
	Real	nInt, xx, yy, aa, bb;
	Real	tIme;

	for (j = 0; j < selectedParticles.size(); j++) {
	// for (j = 0; j < 4; j++) {
		rad0 = pd.getRadius(selectedParticles[j]);
		// rad0 = pd.getRadius(j);
		if ((j % 4) == 0) {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, red);
			// MiniGL::drawSphere(pd.getPosition(j), rad0, red);
		} else if ((j % 4) == 1) {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, grn);
			// MiniGL::drawSphere(pd.getPosition(j), rad0, grn);
		} else if ((j % 4) == 2) {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, blu);
			// MiniGL::drawSphere(pd.getPosition(j), rad0, blu);
		} else if ((j % 4) == 3) {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, blk);
			// MiniGL::drawSphere(pd.getPosition(j), rad0, blk);
		}

		// draw the circles
		/*
		if (j == 0) {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, red);
		} else {
			MiniGL::drawSphere(pd.getPosition(selectedParticles[j]), rad0, grn);
		}
		*/

		// draw the centers (vertices) of the circles
		rad0 = pd.getnInt (j);
		pos = pd.getPosition(selectedParticles[j]);
		MiniGL::drawPoint (pos, 7.5, blu);

		rad0 = pd.getnInt (selectedParticles[j]);
		// rad0 = pd.getnInt (j);
		// fprintf (stderr,"%d %d\n",j,(unsigned int) rad0);
		if (j == 0) {
			for (i = 0 ; i < (unsigned int) rad0 ; i+=2) {
				// pos = pd.getInt (selectedParticles[j], i);
				pos = pd.getInt (j, i);
				// fprintf (stderr,"%d %d / %d %f %f ",j,i,(unsigned int)rad0,pos[0],pos[2]);
				// MiniGL::drawPoint (pos, 7.5, blu);

				// pos3 = pd.getInt (selectedParticles[j], i+1);
				pos3 = pd.getInt (j, i+1);
				// fprintf (stderr,"%f %f\n",pos3[0],pos3[2]);
				// MiniGL::drawPoint (pos3, 7.5, blu);

				// MiniGL::drawVector (pos, pos3, 1.5, blk);
			}
		}
		// pd.setnInt (selectedParticles[j], 0);
	}

	// draw the current network
	if (currentVertices.size() > 0) {
		for (j = 0 ; j < currentVertices.size()-1 ; j++) {
			pos1 = cd.getCurrA(j);
			// pos2 = cd.getCurrB(j);
			// fprintf (stderr,"render: %d %f %f\n",i,pos[0],pos[2]);
			/*
			if (j == 0) {
				MiniGL::drawPoint (pos1, 7.5, grn);
			} else {
				MiniGL::drawPoint (pos1, 7.5, blk);
			}
			*/

			MiniGL::drawPoint (pos1, 7.5, red);
			// MiniGL::drawPoint (pos2, 7.5, blk);
			// MiniGL::drawVector (pos1, pos2, 1.5, red);
		}
	}

	MiniGL::drawTime( TimeManager::getCurrent ()->getTime ());
}


/** Create a particle model mesh 
*/
void createMesh()
{
	TriangleModel::ParticleMesh::UVs uvs;
	uvs.resize(nRows*nCols);

	const Real dy = width / (Real)(nCols - 1);
	const Real dx = height / (Real)(nRows - 1);
	/*
	const Real dy = width / (Real) nCols;
	const Real dx = height / (Real) nRows;
	*/
	Real	rad0;

	Vector3r points[nRows*nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			const Real y = (Real)dy*j;
			const Real x = (Real)dx*i;
			points[i*nCols + j] = Vector3r(x, 1.0, y);

			uvs[i*nCols + j][0] = x/width;
			uvs[i*nCols + j][1] = y/height;
		}
	}
	const int nIndices = 6 * (nRows - 1)*(nCols - 1);

	TriangleModel::ParticleMesh::UVIndices uvIndices;
	uvIndices.resize(nIndices);

	unsigned int indices[nIndices];
	int index = 0;
	for (int i = 0; i < nRows - 1; i++)
	{
		for (int j = 0; j < nCols - 1; j++)
		{
			int helper = 0;
			if (i % 2 == j % 2)
				helper = 1;

			indices[index] = i*nCols + j;
			indices[index + 1] = i*nCols + j + 1;
			indices[index + 2] = (i + 1)*nCols + j + helper;

			uvIndices[index] = i*nCols + j;
			uvIndices[index + 1] = i*nCols + j + 1;
			uvIndices[index + 2] = (i + 1)*nCols + j + helper;
			index += 3;

			indices[index] = (i + 1)*nCols + j + 1;
			indices[index + 1] = (i + 1)*nCols + j;
			indices[index + 2] = i*nCols + j + 1 - helper;

			uvIndices[index] = (i + 1)*nCols + j + 1;
			uvIndices[index + 1] = (i + 1)*nCols + j;
			uvIndices[index + 2] = i*nCols + j + 1 - helper;
			index += 3;
		}
	}
	model.addTriangleModel(nRows*nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);
	
	CurrentData &cd = model.getCurrents();
	ParticleData &pd = model.getParticles();

	// fprintf (stderr,"main: %d %d\n",pd.size(),cd.size());

	// std::cout << "num part " << pd.getNumberOfParticles() << "\n";
	for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
	{
		pd.setMass(i, 1.0f);
		// rad0 = (Real) rand() / RAND_MAX * 2;
		rad0 = (Real) 2;
		// pd.setRadius(i, (Real) rand() / RAND_MAX * 2);
		pd.setRadius(i, rad0);
	}

	// Set mass of points to zero => make it static
	// pd.setMass(0, 20.0);
	/*
	pd.setMass(0, 0.0);
	pd.setMass((nRows-1)*nCols, 0.0);
	*/

	// init constraints
	for (unsigned int cm = 0; cm < model.getTriangleModels().size(); cm++)
	{
		// const unsigned int offset = model.getTriangleModels()[cm]->getIndexOffset();
		// const unsigned int offset = 0.1;
		const unsigned int offset = 0.0;
		IndexedFaceMesh &mesh = model.getTriangleModels()[cm]->getParticleMesh();
		const unsigned int nEdges = mesh.numEdges();
		const IndexedFaceMesh::Edge *edges = mesh.getEdges().data();

		// fprintf (stderr,"init: %d %d\n",cm,nEdges);
		// distance constraints
		for (unsigned int i = 0; i < nEdges; i++)
		{
			const unsigned int v1 = edges[i].m_vert[0] + offset;
			const unsigned int v2 = edges[i].m_vert[1] + offset;

			// fprintf (stderr,"\t%d %d %d\n",i,v1,v2);
			model.addGenericDistanceConstraint(v1, v2);
		}

		// bending constraints
		const unsigned int *tris = mesh.getFaces().data();

		for (unsigned int i = 0; i < nEdges; i++)
		{
			const int tri1 = edges[i].m_face[0];
			const int tri2 = edges[i].m_face[1];
			if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff))
			{
				// Find the triangle points which do not lie on the axis
				const int axisPoint1 = edges[i].m_vert[0];
				const int axisPoint2 = edges[i].m_vert[1];
				int point1 = -1;
				int point2 = -1;
				for (int j = 0; j < 3; j++)
				{
					if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2))
					{
						point1 = tris[3 * tri1 + j];
						break;
					}
				}
				for (int j = 0; j < 3; j++)
				{
					if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2))
					{
						point2 = tris[3 * tri2 + j];
						break;
					}
				}
				if ((point1 != -1) && (point2 != -1))
				{
					const unsigned int vertex1 = point1 + offset;
					const unsigned int vertex2 = point2 + offset;
					const unsigned int vertex3 = edges[i].m_vert[0] + offset;
					const unsigned int vertex4 = edges[i].m_vert[1] + offset;
					model.addGenericIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
				}
			}
		}
	}

	std::cout << "Number of triangles: " << nIndices / 3 << "\n";
	std::cout << "Number of vertices: " << nRows*nCols << "\n";

}

void TW_CALL setTimeStep(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	TimeManager::getCurrent()->setTimeStepSize(val);
}

void TW_CALL getTimeStep(void *value, void *clientData)
{
	*(Real *)(value) = TimeManager::getCurrent()->getTimeStepSize();
}

void TW_CALL setVelocityUpdateMethod(const void *value, void *clientData)
{
	const short val = *(const short *)(value);
	((TimeStepController*)clientData)->setVelocityUpdateMethod((unsigned int)val);
}

void TW_CALL getVelocityUpdateMethod(void *value, void *clientData)
{
	*(short *)(value) = (short)((TimeStepController*)clientData)->getVelocityUpdateMethod();
}

void TW_CALL setStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((GenericConstraintsModel*)clientData)->setClothStiffness(val);
}

void TW_CALL getStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((GenericConstraintsModel*)clientData)->getClothStiffness();
}

void TW_CALL setBendingStiffness(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	((GenericConstraintsModel*)clientData)->setClothBendingStiffness(val);
}

void TW_CALL getBendingStiffness(void *value, void *clientData)
{
	*(Real *)(value) = ((GenericConstraintsModel*)clientData)->getClothBendingStiffness();
}
