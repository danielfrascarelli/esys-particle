#!/bin/sh
# A script that runs all the scripts to test ESys-Particle installation

START=`date +%s`
 
python out2file.py w $1"-Test.txt" "# Testing build: "$1

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/bingle.py
then python out2file.py a $1"-Test.txt" "bingle.py pass"
else
	python out2file.py a $1"-Test.txt" "bingle.py FAIL"
fi

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/bingle_output.py
then python out2file.py a $1"-Test.txt" "bingle_output.py pass"
else
	python out2file.py a $1"-Test.txt" "bingle_output.py FAIL"
fi

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/bingle_chk.py
then python out2file.py a $1"-Test.txt" "bingle_chk.py pass"
else
	python out2file.py a $1"-Test.txt" "bingle_chk.py FAIL"
fi
sh move_output_files.sh bingle_chk bingle_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/bingle_vis.py
then python out2file.py a $1"-Test.txt" "bingle_vis.py pass"
else
	python out2file.py a $1"-Test.txt" "bingle_vis.py FAIL"
fi
sh move_output_files.sh bingle_vis bingle_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/gravity.py
then python out2file.py a $1"-Test.txt" "gravity.py pass"
else
	python out2file.py a $1"-Test.txt" "gravity.py FAIL"
fi
sh move_output_files.sh gravity gravity_data

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/gravity_cube.py
then python out2file.py a $1"-Test.txt" "gravity_cube.py pass"
else
	python out2file.py a $1"-Test.txt" "gravity_cube.py FAIL"
fi
sh move_output_files.sh gravity_cube gravity_data

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/slope_fail.py
then python out2file.py a $1"-Test.txt" "slope_fail.py pass"
else
	python out2file.py a $1"-Test.txt" "slope_fail.py FAIL"
fi
sh move_output_files.sh slope_fail slope_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/slope_friction.py
then python out2file.py a $1"-Test.txt" "slope_friction.py pass"
else
	python out2file.py a $1"-Test.txt" "slope_friction.py FAIL"
fi
sh move_output_files.sh slope_friction slope_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/slope_friction_floor.py
then python out2file.py a $1"-Test.txt" "slope_friction_floor.py pass"
else
	python out2file.py a $1"-Test.txt" "slope_friction_floor.py FAIL"
fi
sh move_output_files.sh slope_friction_floor slope_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/slope_friction_walls.py
then python out2file.py a $1"-Test.txt" "slope_friction_walls.py pass"
else
	python out2file.py a $1"-Test.txt" "slope_friction_walls.py FAIL"
fi
sh move_output_files.sh slope_friction_walls slope_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/hopper_flow.py
then python out2file.py a $1"-Test.txt" "hopper_flow.py pass"
else
	python out2file.py a $1"-Test.txt" "hopper_flow.py FAIL"
fi
sh move_output_files.sh hopper_flow flow_data_

if mpirun -np 2 '/home/jrahardjo/BUILD/esys-particle/install/bin/esysparticle' Scripts/rot_compress.py
then python out2file.py a $1"-Test.txt" "rot_compress.py pass"
else
	python out2file.py a $1"-Test.txt" "rot_compress.py FAIL"
fi
mkdir -p rot_compress
mv *.dat rot_compress

END=`date +%s`
ELAPSED=$(( $END - $START ))

python out2file.py a $1"-Test.txt" "Time consumed by testing: "$ELAPSED" seconds"