#include <gtest/gtest.h>
#include <ComanODE/ComanSimClass.h>

ComanSimClass odesim;
dsFunctions fn;

// #include <unistd.h>
// #include <sys/types.h>
// #include <pwd.h>
// const char *homedir;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
	odesim.nearCallback(data, o1, o2);
}

void simLoop(int pause)
{
	odesim.simLoop(pause);

// #define TEST_ODE_TIME
// #ifdef TEST_ODE_TIME
// 	static int count = 1;
// 	static double aveT = 0.0;
// 	double t1 = get_time();
// #endif

	dSpaceCollide(odesim.getSpaceID(), 0, &nearCallback);
	odesim.UpdateODE();

// #ifdef TEST_ODE_TIME
// 	double t2 = get_time();
// 	aveT += (t2 - t1);
// 	COUT("ODE simulation step using ", t2 - t1, " s. Average time is ", aveT / count);
// 	count ++;
// #endif
}

void command(int cmd)
{
	odesim.command(cmd);
}

void start()
{
	odesim.start();
}

void setDrawStuff(dsFunctions &fn)
{
	// auto f1 = boost::bind(&ODESimBaseClass::start, this);
	fn.version = DS_VERSION;
	fn.start   = &start;
	fn.step    = &simLoop;
	fn.command = &command;
#ifdef USE_WINDOWS
	fn.path_to_textures = ".";
#else
	// if ((homedir = getenv("HOME")) == NULL) {
	// 	homedir = getpwuid(getuid())->pw_dir;
	// }

	const char * pathname = ROBOT_PARA_YAML_PATH;
	const char * pathname_postfix = "../../..";
	char * textures_path = (char *) malloc(1 + strlen(pathname) + strlen(pathname) );
	strcpy(textures_path, pathname);
	strcat(textures_path, pathname_postfix);
	fn.path_to_textures = textures_path;

	// fn.path_to_textures = "/home/czhou/src/ode-0.14/drawstuff/textures";
#endif
}

TEST(unittest, ODEtest)
{
	setDrawStuff(fn);
	EXPECT_EQ(0, odesim.main(fn));
}

// TEST(unittest, SVM_collect_train_data)
// {
// 	setDrawStuff(fn);

// 	RobotParaClass::robot_push_force[1] = 0;
// 	// x fall
// 	std::vector<double> v1 = {0, 10, 20, 20, 25, 25, 30, 30};
// 	double push = 130;

// 	for (int i = 0; i < v1.size(); ++i) {
// 		push += v1[i];
// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[0] = push;
// 		ASSERT_EQ(0, odesim.main(fn));

// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[0] = -push;
// 		ASSERT_EQ(0, odesim.main(fn));
// 	}

// 	// x stable
// 	std::vector<double> v2 = {0, 10, 20, 20, 30};
// 	push = 130;

// 	for (int i = 0; i < v2.size(); ++i) {
// 		push -= v2[i];
// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[0] = push;
// 		ASSERT_EQ(0, odesim.main(fn));

// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[0] = -push;
// 		ASSERT_EQ(0, odesim.main(fn));
// 	}

// 	RobotParaClass::robot_push_force[0] = 0;
// 	// y fall
// 	std::vector<double> v3 = {0, 10, 20, 20, 25, 25, 30, 30};
// 	push = 200;

// 	for (int i = 0; i < v3.size(); ++i) {
// 		push += v3[i];
// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[1] = push;
// 		ASSERT_EQ(0, odesim.main(fn));

// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[1] = -push;
// 		ASSERT_EQ(0, odesim.main(fn));
// 	}

// 	// y stable
// 	std::vector<double> v4 = {0, 10, 20, 20, 30, 30};
// 	push = 200;

// 	for (int i = 0; i < v4.size(); ++i) {
// 		push -= v4[i];
// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[1] = push;
// 		ASSERT_EQ(0, odesim.main(fn));

// 		odesim = ComanSimClass();
// 		RobotParaClass::robot_push_force[1] = -push;
// 		ASSERT_EQ(0, odesim.main(fn));
// 	}

// }

int main(int argc, char* argv[])
{
	// if folder does not exist, mkdir
	struct stat st = {0};
	if (stat("/tmp/mrdplot", &st) == -1) {
		mkdir("/tmp/mrdplot", 0700);
	}
	// testing::AddGlobalTestEnvironment(new FooEnvironment); // 添加全局事件
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}