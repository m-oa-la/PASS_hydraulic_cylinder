'use strict';

////////////////////////////////////////////////////////////////
////////////////////////// VARIABLES ///////////////////////////
////////////////////////////////////////////////////////////////

var icon = 'imgs/icon-hcl.png';
var menuIcon = 'imgs/icon-hcl.png';
var flanged_endEye = 'imgs/flanged.png';
var flanged_located_at_end_cover_endEye = 
  'imgs/flanged_located_at_cylinder_tube.png';
var flanged_located_at_end_cover_thread = 
  'imgs/flanged_located_at_cylinder_tube_thread.png';
var flanged_rod_always_inside_cylinder = 
  'imgs/flanged_rod_always_inside_cylinder.png';
var flanged_rod_always_inside_cylinder_thread = 
  'imgs/flanged_rod_always_inside_cylinder_thread.png';
var endEye_endEye = 'imgs/measures_2.png';
var trunnion_endEye = 'imgs/trunnion_endeye.png';
var trunnion_thread = 'imgs/trunnion_thread.png';
var flanged_thread = 'imgs/flange_threded.png';
var endEye_thread = 'imgs/endeye_threaded.png';
var trunnion_weld_partial = 'imgs/partial_trunnion.png';
var trunnion_weld_fillet = 'imgs/fillettrunnion.png';
var steering_edge = 'imgs/steeringedge.png';
var flange_fillet = 'imgs/fillet_flange.png';
var flange_partial = 'imgs/flanged_partial.png';
var flange_partial_mono = 'imgs/flange_partial_mono.png';
var flanged_end = 'imgs/flanged_end2.png';
var ring = 'imgs/ring.png';
var du = 'imgs/du.png';

var MatList = ["S355J2+N, EN10025-3", 
	"E355+SR, EN10305-1", 
	"1.4462, EN 10088-3", 
	"1.4460, EN 10088-3", 
	"1.4418 +QT760 EN 10088-3", 
	"1.4418 +QT900 EN 10088-3", 
	"1.4404 EN 10088-3", 
	"1.4403 EN 10088-3"];

var f_y_S_355 = [[0, 16, 355, 470], 
	[16, 40, 345, 470], 
	[40, 63, 335, 470], 
	[63, 80, 325, 470], 
	[80, 100, 315, 470], 
	[100, 150, 295, 450], 
	[150, 200, 285, 450], 
	[200, 250, 285, 450], 
	[250, 400, 265, 450]];

var f_y_E_355 = [[0, 160, 450, 580], 
	[160, 400, 420, 580]];

var f_y_1_4462 = [[0, 16, 650, 850], 
	[16, 160, 450, 650]];

var f_y_1_4460 = [[0, 10, 610, 770], 
	[10, 16, 560, 770], 
	[16, 160, 460, 620]];

var f_y_1_4418QT760 = [[0, 400, 550, 760]];

var f_y_1_4418QT900 = [[0, 400, 700, 900]];

var f_y_1_4404 = [[0, 400, 235, 500]];

var f_y_1_4403 = [[0, 400, 225, 500]];

////////////////////////////////////////////////////////////////
////////////////////// UTILITY FUNCTIONS ///////////////////////
////////////////////////////////////////////////////////////////

function roundToDecimal(number, precision) {
  /* Rounds a number to the closest decimal number with a given precision.
  :arg number: The number to be rounded.
  :arg precision: The number of decimals to be included in the rounding.
  :return: A float. */
  precision = Math.round(precision);

  if (precision < 0){
    return number;
  }
  else
  {
    return Math.round(Math.pow(10, precision) * number) / Math.pow(10, precision);
  }
} // End of roundToDecimal function

function find_minimum(arr) {
	if (arr.length == 0) {
		return null;
	}

	var min = Infinity;
	var x = undefined;
	for (var i = 0; i < arr.length; i++) {
		x = arr[i];
		if (x < min) {
			min = x;
		}
	}
	return min;
}

function calcArea(outerD, innerD) {
	/* Calculates the area of a circular cross section. 
	 * arg outerD: outer diameter
	 * arg innerD: inner diameter */
	var pi = Math.PI;
	return pi/4*(Math.pow(outerD,2) - Math.pow(innerD,2));
}

function calcInertiaMoment(outerDiameter, innerDiameter) {
  /* Calculates the second moment of inertia for a circular cross section.
  :arg outer_diameter: The outer diameter of the cross section.
  :arg inner_diameter: The inner diameter of the cross section. */
  return (Math.PI/64) * (Math.pow(outerDiameter,4) - Math.pow(innerDiameter,4));
} // End of calcInertiaMoment function

function calcZFactor(I_tube, I_rod, L_euler, L1_Le, L2max_Le) {
  /* Calculates the Z-factor used in buckling calculations.
  :arg I_tube: The second moment of inertia of the tube.
  :arg I_rod: The second moment of inertia of the rod.
  :arg L_euler: The Euler length of the rod.
  :arg L1_Le: The ratio between cylinder tube and the Euler length.
  :arg L2max_Le: The ratio between the fully extract rod length and the Euler length. */
  var Z = L1_Le / I_tube 
    + L2max_Le / I_rod 
    + (1 / I_tube - 1 / I_rod) 
    * (L_euler / (2*Math.PI)) 
    * Math.sin(2 * Math.PI * L1_Le / L_euler);

  return Z;
} // End of calcZFactor

function calc_C2(DI, es, v, fmin, P) {
  /* Calculates the C2 factor used in the end cover thickness rules.
  :arg DI: Inner diameter of the cylinder tube.
  :arg es: End cover thickness (?).
  :arg v: Poisson ratio.
  :arg fmin: Minimum stress.
  :arg P: Pressure in MPa.
  :return: A float. */
  P = P / 10; // Calculations assume bar

  var g = DI/(DI+es);
  var H = Math.pow(12*(1 - v*v), 0.25) * Math.pow(es/(DI+es), 0.5);
  var J = 3*fmin/P - DI*DI/(4*es*(DI+es)) - 1;
  var U = 2*(2 - v*g)/Math.sqrt(3*(1 - v*v));
  var f1 = 2*g*g - g*g*g*g;

  var A = (3*U*DI/(4*es) - 2*J)*(1+v)*(1 + (1-v)*es/(DI + es));
  var B = ((3*U*DI/(8*es) - J)*H*H - 3/2*(2-v*g)*g)*H;
  var F = (3/8*U*g + 3/16*f1*(DI+es)/es - 2*J*es/(DI+es))*H*H - 3*(2-v*g)*g*es/(DI+es);
  var G = (3 / 8 * f1 - 2 * J * Math.pow(es / (DI + es), 2)) * H;

  var a = B/A;
  var b = F/A;
  var c = G/A;

  var N = b/3 - a*a/9;
  var Q = c/2 - a*b/6 + a*a*a/27;
  var K = N*N*N/(Q*Q);
  var S = 0;

  if (Q >= 0) {
    S = Math.pow(Q*(1 + Math.sqrt(1+K)), 1/3);
  } else {
    S = -Math.pow(Math.abs(Q)*(1 + Math.sqrt(1+K)), 1/3);
  }

  var C2 = (DI+es)*(N/S - S - a/3)/(DI*Math.sqrt(P/fmin));
  return Math.abs(C2);
} // End of calc_C2 function

function getYieldAndStrength(thickness, material) {
	/* Finds the yield stress and tensile strength of an circular 
	cross section based on the material type, shell thickness and 
	material.
	:arg thickness: The shell thickness. (OD - ID) / 2.
	:arg material: The material properties.
	:return: list of three floats. */
	var def = material.def;
  	var yieldTable = material.f_y;
	var f_y = undefined;
	var f_yt = undefined;
	var TS = undefined;
	if (def == "standard") {
		f_y = yieldTable[0][2];
		f_yt = Infinity;
		TS = yieldTable[0][3];
		for (var j = 0; j < yieldTable.length; j++) {
			if (thickness >= yieldTable[j][0] && thickness < yieldTable[j][1]) {
				f_y = yieldTable[j][2];
				TS = yieldTable[j][3];
			}
		}
	} else if (def == "custom") {
		f_y = yieldTable[0][2];
		f_yt = yieldTable[0][3];
		TS = yieldTable[0][4];
		for (var j = 0; j < yieldTable.length; j++) {
			if (thickness >= yieldTable[j][0] && thickness < yieldTable[j][1]) {
				f_y = yieldTable[j][2];
				f_yt = yieldTable[j][3];
				TS = yieldTable[j][4];
			}
		}
	}
  	return [f_y, f_yt, TS];
} // End of getYieldAndStrength function

function calc_nominal_stress(thickness, factors, material) {
	/* Calculates the nominal design stress based on the material yield stress
	 * and tensile strength at room temperature and the material yield stress at
	 * design temperature.
	:arg thickness: The shell thickness.
	:arg factors: Nominal stress weighting factors.
	:arg material: The material properties. */
	if (factors.length != 3) {
		return;
	}

	var data = getYieldAndStrength(thickness, material)
	var fracs = data.map((e,i) => e / factors[i]);
	var sigma = find_minimum(fracs);
	return sigma;
} // End of calc_nominal_stress function

function configureCylinder(handle) {
	// Configure the cylinder parts
	var cylinder = {};
	cylinder["tube"] = configureTube(handle);
	cylinder["rod"] = configureRod(handle);
	cylinder["endCover"] = configureEndCover(handle);
	cylinder["piston"] = configurePiston(handle);
	cylinder["stuffingBox"] = configureStuffingBox(handle);

	// Configure non-part specific information
	cylinder["maxLength"] = cylinder.tube.length + cylinder.rod.maxLength;
	cylinder["minLength"] = cylinder.tube.length + cylinder.rod.minLength;
	cylinder["euler"] = handle.get("eL");
	cylinder["mass"] = handle.get("m_cyl");
	cylinder["corrisonAllowance"] = 0.3;
	cylinder["medium"] = handle.get("medium_type");
	cylinder["application"] = handle.get("pv_func");
	
	// Calculation of cylinder force
	cylinder["pullPressure"] = handle.get("DPpull");
	cylinder["pushPressure"] = handle.get("DPpush");
	var pullPressure = cylinder.pullPressure;
	var pushPressure = cylinder.pushPressure;
	var area = cylinder.piston.area * Math.pow(10, -6); // Convert to SI
	var pressure = 0;
	if (cylinder.tube.end.typeForce == "trunnion" && cylinder.rod.inside == "yes") {
		pressure = Math.max(pullPressure, pushPressure) * Math.pow(10, 5); // Convert to SI
	} else {
		pressure = pushPressure * Math.pow(10, 5); // Convert to SI
	}
	cylinder["force"] = pressure * area;
	
	// Buckling safety factor calculation method
	cylinder["bucklingMethod"] = handle.get("calcMethod");
	if (typeof cylinder.bucklingMethod != "undefined") {
		if (cylinder.bucklingMethod.indexOf("force") != -1) {
			cylinder["forceCurve"] = handle.get("fCurve");
		} else if (cylinder.bucklingMethod.indexOf("pressure") != -1) {
			cylinder["pressureCurve"] = handle.get("pCurve");
		}
	}
	return cylinder;
}

function configureTube(handle) {
	/* Configures a tube object.
	 * arg handle: object, page handle
	 * return: object */
	var tube = {}
	tube["material"] = handle.get("tubeMaterial");
	tube["length"] = handle.get("L1");
	tube["outerD"] = handle.get("OD");
	tube["innerD"] = handle.get("DI");
	tube["thickness"] = (tube.outerD - tube.innerD)/2;
	var materialData = getYieldAndStrength(tube.thickness, tube.material);
	tube["yield"] = materialData[0];
	tube["yieldDesign"] = materialData[1];
	tube["tensile"] = materialData[2];
	tube["weldEfficiency"] = handle.get("v");
	tube["tolerance"] = handle.get("cm_shell");
	tube["euler"] = handle.get("eL");
	tube["area"] = calcArea(tube.outerD, tube.innerD);
	tube["inertia"] = calcInertiaMoment(tube.outerD, tube.innerD);
	tube["end"] = configureTubeEnd(handle);
	return tube;
}

function configureTubeEnd(handle) {
	/* Configures a tube end object.
	 * arg handle: object, page handle
	 * return: object */
	var tubeEnd = {};
	tubeEnd["type"] = handle.get("tubeEnd"); // endEye, trunnion, flanged, flanged_located_at_end_cover
	if (tubeEnd["type"] == "endEye") {
		tubeEnd["width"] = handle.get("T_tube");
		tubeEnd["outerR"] = handle.get("R_tube");
		tubeEnd["innerD"] = handle.get("d_tube");
		tubeEnd["material"] = handle.get("tubeEndEyeMaterial");
		tubeEnd["attachment"] = handle.get("tubeEndEyeAttachment");
		if (tubeEnd["attachment"] == "machined") {
			
		} else if (tubeEnd["attachment"] == "welded") {
			tubeEnd["weldType"] = handle.get("tubeEndEyeWeldType");
			if (tubeEnd["weldType"] == "pPen" || tubeEnd["weldType"] == "fWeld") {
				tubeEnd["fatigue"] = handle.get("fatigueChoice1");
				if (tubeEnd["fatigue"] == "fChoice2") {
					tubeEnd["weldArea"] = handle.get("weld_area_tube");
					tubeEnd["fatigueCycles"] = handle.get("n_cycles_manufacturer");
				}
			}
		}
	} else if (tubeEnd["type"] == "trunnion") {
		tubeEnd["weldType"] = handle.get("trunnionWeldType");
		if (tubeEnd["weldType"] == "pPen") {
			tubeEnd["fatigue"] = handle.get("fatigueChoice3");
			if (tubeEnd["fatigue"] == "fChoice2") {
				tubeEnd["weldMinLeg"] = handle.get("weld_throat_trunnion");
				tubeEnd["fatigueCycles"] = handle.get("n_cycles_manufacturer_trunnion");
			}
		} else if (tubeEnd["weldType"] == "fWeld") {
			tubeEnd["fatigue"] = handle.get("fatigueChoice3");
			if (tubeEnd["fatigue"] == "fChoice2") {
				tubeEnd["weldThroat"] = handle.get("weld_throat_trunnion");
				tubeEnd["fatigueCycles"] = handle.get("n_cycles_manufacturer_trunnion");
			}
		}
	} else if (tubeEnd["type"] == "flanged" || tubeEnd["type"] == "flanged_located_at_end_cover") {
		tubeEnd["flangeThickness"] = handle.get("T_flange");
		tubeEnd["boltN"] = handle.get("n_bolt_flange");
		tubeEnd["boltD"] = handle.get("D_bolt_flange");
		tubeEnd["material"] = handle.get("flangeMaterial");
		tubeEnd["weldType"] = handle.get("flangedWeldType");
		if (tubeEnd["weldType"] == "pPen" || tubeEnd["weldType"] == "pPen1") {
			tubeEnd["fatigue"] = handle.get("fatigueChoice4");
			tubeEnd["weldMinLeg"] = handle.get("weld_leg_flanged");
			if (tubeEnd["fatigue"] == "fChoice2") {
				tubeEnd["fatigueCycles"] = handle.get("n_cycles_manufacturer_flanged");
			}
		} else if (tubeEnd["weldType"] == "fWeld") {
			tubeEnd["fatigue"] = handle.get("fatigueChoice4");
			tubeEnd["weldThroat"] = handle.get("weld_throat_flanged");
			if (tubeEnd["fatigue"] == "fChoice2") {
				tubeEnd["fatigueCycles"] = handle.get("n_cycles_manufacturer_flanged");
			}
		}
	}

	var calcMethod = handle.get("calcMethod");
	if (calcMethod == "accurate_static" || calcMethod == "accurate_force" || calcMethod == "accurate_pressure"){
		tubeEnd["friction"] = handle.get("my_end_eye");
	}
	return tubeEnd;
}

function configureRod(handle) {
	/* Configures a rod object.
	 * arg handle: object, page handle
	 * return: object */
	var rod = {};
	rod["material"] = handle.get("rodMaterial");
	rod["maxLength"] = handle.get("L2max");
	rod["minLength"] = rod.maxLength - handle.get("Stroke");
	rod["outerD"] = handle.get("OD_rod");
	rod["innerD"] = handle.get("DI_rod");
	rod["thickness"] = (rod.outerD - rod.innerD)/2;
	var materialData = getYieldAndStrength(rod.thickness, rod.material);
	rod["yield"] = materialData[0];
	rod["yieldDesign"] = materialData[1];
	rod["tensile"] = materialData[2];
	rod["euler"] = handle.get("eL");
	rod["area"] = calcArea(rod.outerD, rod.innerD);
	rod["inertia"] = calcInertiaMoment(rod.outerD, rod.innerD);
	rod["end"] = configureRodEnd(handle);

	if (handle.get("tubeEnd") == "trunnion") {
		rod["inside"] = handle.get("rodInside");
	}

	return rod;
}

function configureRodEnd(handle) {
	/* Configures a rod end object.
	 * arg handle: object, page handle
	 * return: object */
	var rodEnd = {};
	rodEnd["type"] = handle.get("rodEnd");
	
	if (rodEnd["type"] == "endEye") {
		rodEnd["material"] = handle.get("rodEndEyeMaterial");
		rodEnd["width"] = handle.get("T_rod");
		rodEnd["outerR"] = handle.get("R_rod");
		rodEnd["innerD"] = handle.get("d_rod");
		rodEnd["attachment"] = handle.get("rodEndEyeAttachment");
		if (rodEnd["attachment"] == "threaded") {
			rodEnd["eyeThreadD"] = handle.get("Md_rod_eye");
			rodEnd["eyeThreadP"] = handle.get("xP_rod_eye");
			rodEnd["eyeThreadL"] = handle.get("Le_rod_eye");
		} else if (rodEnd["attachment"] == "welded") {
			rodEnd["weldType"] = handle.get("rodEndEyeWeldType");
			if (rodEnd["weldType"] == "pPen" || rodEnd["weldType"] == "fWeld") {
				rodEnd["fatigue"] = handle.get("fatigueChoice2");
				if (rodEnd["fatigue"] == "fChoice2") {
					rodEnd["weldArea"] = handle.get("weld_area_rod");
					rodEnd["fatigueCycles"]= handle.get("n_cycles_manufacturer_rod_eye");
				}
			}
		}
	} else if (rodEnd["type"] == "threaded_fixed" || rodEnd["type"] == "threaded_pinned") {
		rodEnd["threadD"] = handle.get("Md_rod");
		rodEnd["threadP"] = handle.get("xP_rod");
		rodEnd["threadL"] = handle.get("Le_rod");
	}
	
	var calcMethod = handle.get("calcMethod");
	if (calcMethod == "accurate_static" || calcMethod == "accurate_force" || calcMethod == "accurate_pressure"){
		rodEnd["friction"] = handle.get("my_end_eye");
	}

	return rodEnd;
}

function configureEndCover(handle) {
	/* Configures a end cover object.
	 * arg handle: object, page handle
	 * return: object */
	var endCover = {};
	endCover["material"] = handle.get("ecMaterial");
	endCover["thickness"] = handle.get("ect");
	endCover["attachment"] = handle.get("ecAttachment");
	if (endCover.attachment == "welded") {
		endCover["steeringLength"] = handle.get("L_steering_edge");
	} else if (endCover.attachment == "bolted") {
		endCover["bolt"] = configureEndCoverBolt(handle);
	}
	return endCover;
}

function configureEndCoverBolt(handle) {
	/* Configures a end cover bolt object.
	 * arg handle: object, page handle
	 * return: object */
	var bolt = {};
	bolt["material"] = handle.get("ecBoltMaterial");
	bolt["nominalD"] = handle.get("d_n_bolt_ec");
	bolt["number"] = handle.get("n_bolts_ec");
	bolt["threadP"] = handle.get("L_p_bolt_ec");
	bolt["threadL"] = handle.get("l_thread_bolt_ec");
	bolt["stressA"] = handle.get("A_s_bolt_ec");
	bolt["circleD"] = handle.get("C_bolt_circle");
	return bolt;
}

function configurePiston(handle) {
	/* Configures a piston object.
	 * arg handle: object, page handle
	 * return: object */
	var piston = {};
	piston["material"] = handle.get("pistonMaterial");
	piston["stroke"] = handle.get("Stroke");
	piston["guidingL"] = handle.get("L3");
	piston["threadD"] = handle.get("Md_piston");
	piston["threadP"] = handle.get("xP_piston");
	piston["threadL"] = handle.get("Le_piston");
	var outerD = handle.get("DI");
	var innerD = 0;
	if (handle.get("tubeEnd") == "trunnion") {
		if (handle.get("rodInside")) {
			innerD = handle.get("DI_rod");
		}
	} 
	piston["outerD"] = outerD;
	piston["innerD"] = innerD;
	piston["area"] = calcArea(outerD, innerD);
	return piston;
}

function configureStuffingBox(handle) {
	/* Configures a stuffing box object.
	 * arg handle: object, page handle
	 * return: object */
	var stuffingBox = {};
	stuffingBox["material"] = handle.get("sbMaterial");
	stuffingBox["attachment"] = handle.get("sbAttachment");
	stuffingBox["threadType"] = handle.get("sbThreadType");
	if (stuffingBox.attachment == "threaded") {
		stuffingBox["threadMetricD"] = handle.get("Md_sb");
		stuffingBox["threadP"] = handle.get("xP_sb");
		stuffingBox["threadL"] = handle.get("Le_sb");
		stuffingBox["threadReliefD"] = handle.get("Du");
	} else if (stuffingBox.attachment == "bolted") {
		stuffingBox["bolt"] = configureStuffingBoxBolt(handle);
	}
	return stuffingBox;
}

function configureStuffingBoxBolt(handle) {
	/* Configures a stuffing box bolt object.
	 * arg handle: object, page handle
	 * return: object */
	var bolt = {};
	bolt["material"] = handle.get("sbBoltMaterial");
	bolt["number"] = handle.get("n_bolts");
	bolt["nominalD"] = handle.get("d_n_bolt");
	bolt["threadP"] = handle.get("L_p_bolt");
	bolt["stressA"] = handle.get("A_s_bolt");
	bolt["threadL"] = handle.get("l_thread_bolt")
	return bolt;
}

////////////////////////////////////////////////////////////////
//////////////////////// RULE FUNCTIONS ////////////////////////
////////////////////////////////////////////////////////////////

function rule_thread_stress(set, error, warn, part_name, Md, xP, Le, DI, OD, material1, 
	material2, P, red, isPiston, isPull) {
	/* Takes in the parameters of a thread and calculates the shear and compression
	stress in it. The function shows an error if the stresses are not according
	to the rules and a warning if they are.
	:arg set: Pointer for setting variables.
	:arg error: Pointer for setting errors.
	:arg warn: Pointer for setting warnings.
	:arg part_name: The name of the part.
	:arg Md: Thread metric diameter.
	:arg xP: Thread pitch.
	:arg Le: Thread engagement length.
	:arg DI: The inner diameter of the cylinder tube.
	:arg OD: The outer diameter of the cylinder tube.
	:arg material1: Material of part 1.
	:arg material2: Material of part 2.
	:arg P: Design pressure, maximum of design pull or push pressure.
	:arg red: Reduction factor (?).
	:arg isPiston: Whether the part is the piston or not.
	:arg isPull: Whether the thread is designed for ONLY pulling or not. */
	var D2 = Md - 0.6495 * xP; // basic pitch diameter of internal thread ISO 724
	var d2 = Md - 0.6495 * xP; // basic pitch diameter of external thread ISO 724
	var D1 = Md - 1.0825 * xP; // basic minor diameter of internal thread ISO 724
	var d1 = Md - 1.0825 * xP; // basic minor diameter of external thread ISO 724

	var fsh = P*(DI*DI - isPull*(d2*d2*red + (1 - red)*OD*OD))/(2*Le*D2)/10;
	var fcomp = P * xP / (Le * (Md * Md - D1 * D1)) / 10
		* (DI * DI - isPull * (d2 * d2 * (isPiston == true ? 1 : 0) 
		+ (isPiston == true ? 0 : 1) * OD * OD)); // unit: N/
	set(part_name + "_fsh", fsh);
	set(part_name + "_comp", fcomp);

	// TODO: Check if yield at room temperature should be added here!
	var mat = Math.min(material1[0] / 1.5, material1[2] / 2.4, 
		material2[0] / 1.5, material2[2] / 2.4);
	var max_shear = 0.8 * mat;
	var mat2 = Math.min(material1[0], material2[0]);

	fsh = roundToDecimal(fsh, 3);
	fcomp = roundToDecimal(fcomp, 3);
	max_shear = roundToDecimal(max_shear, 3);
	mat2 = roundToDecimal(mat2, 3);

	set(part_name + "_f", mat);

	if (!(fsh < max_shear)) {
		error([], "According to $ref the shear stress in the " + part_name 
		+ " thread is too high: \n" + fsh + " MPa > " + max_shear + " MPa.");
	} else {
		warn([], "According to $ref the shear stress in the " + part_name 
		+ " thread is lower than or equal the maximal limit: \n" + fsh + " MPa =< " 
		+ max_shear + " MPa.");
	}
	if (!(fcomp < mat2)) {
		error([], "According to $ref the compressive stress in the " + part_name 
		+ " thread is too high: \n" + fcomp + " MPa > " + mat2 + " MPa.");
	} else {
		warn([], "According to $ref the compressive stress in the " + part_name 
		+ " thread is lower than or equal the maximal limit: \n" + fcomp + " MPa =< " 
		+ mat2 + " MPa.");
	}
} // End of rule_thread_stress function

// FATIGUE RULES
function calc_fatigue(set, error, warn, part_name, F_push, F_pull, weld_type,
  fatigue_type, weld_area, n_cycles_manufacturer) {
	/* Calculates the predicted number of fatigue cycles before failure and compares it
	 * to the number of fatigue cycles before failure provided by the manufacturer.
	 * Displays an error of the rule condition is not met and a warning otherwise.
	 * arg set: Pointer for setting variables.
	 * arg error: Pointer for setting errors.
	 * arg warn: Pointer for setting warnings.
	 * arg part_name: The name of the part.
	 * arg F_push: The push force in kN.
	 * arg F_pull: The pull force in MN.
	 * arg weld_type: The weld type.
	 * arg fatigue_type: The type of fatigue analysis.
	 * arg weld_area: The area of the weld.
	 * arg n_cycles_manufacturer: The number of fatigue cycles provided by manufacturer. */
	if (weld_type != "fPen") {
		if (fatigue_type == "fChoice2") {
			var sigma_weld_pull = F_pull / weld_area;
			var sigma_weld_push = 1000 * F_push / weld_area;

			// Both Fillet weld and partial penetration weld are both considered as group W3 
			// according to DNVGL-RP-C203
			var delta_stress = sigma_weld_pull + sigma_weld_push;
			set("delta_stress_" + part_name, delta_stress);
			var m2 = 3;

			var log_a = 10.97;
			var log_N = 10.97 - 3 * Math.log10(delta_stress);
			var n_fatigue = Math.pow(10, log_N);

			set("sigma_weld_pull_" + part_name, sigma_weld_pull);
			set("sigma_weld_push_" + part_name, sigma_weld_push);
			set("n_cycles_" + part_name, n_fatigue);
			n_fatigue= Math.round(n_fatigue);

			if (n_fatigue < n_cycles_manufacturer)
			{
				error([], "The " + part_name + " weld has a predicted number of cycles to failure of: " 
				  + n_fatigue + " < " + n_cycles_manufacturer + ". Please refer to DNVGL-RP-C203.");
			}
			else
			{
				warn([], "The " + part_name + " weld has a predicted number of cycles to failure of: " 
				  + n_fatigue + " >= " + n_cycles_manufacturer + ". Please refer to DNVGL-RP-C203.");
			}
		} else {
			warn([], "Limit cylinder to 15 000 cycles over life span. Please refer to DNVGL-RP-C203.");
		}
	}
} // End of fatigue function

// BOLT RULES
function calc_bolts(set, error, warn, part_name, d_n_bolt, L_p_bolt, 
  bolt_material, cylinder_material, A_s_bolt, n_bolts, l_thread_bolt, attachment, F_pull, 
  tube_end, Pa) {
	/* Calculates the required engagement length and the maximum shear stress in the 
	bolt. If the requirements are not met an error will be displayed, else a warning
	will be shown.
	:arg set: Pointer for setting variables.
	:arg error: Pointer for setting errors.
	:arg warn: Pointer for setting warnings.
	:arg part_name: The name of the part.
	:arg d_n_bolt: The nominal diameter of the bolts.
	:arg L_p_bolt: The thread pitch of the bolts.
	:arg bolt_material: The bolt material.
	:arg cylinder_material: The cylinder material.
	:arg A_s_bolt: The shear area of the bolts.
	:arg n_bolts: The number of bolts.
	:arg l_thread_bolt: The engagement length of the bolts.
	:arg attachment: The attachment type. Bolted or threaded
	:arg F_pull: The pull force on the part in MN.
	:arg tube_end: The type of tube end.
	:arg Pa: The push force on the part in MN. */
	if (attachment != "bolted") {
		return;
	}

	if (bolt_material.mat == "8.8") {
		var Rp1 = 640;
		var Rp1t = 640;
		var Rm = 800;
	} else if (bolt_material.mat == "10.9") {
		var Rp1 = 900;
		var Rp1t = 900;
		var Rm = 1040;
	}

	var l_thread_bolt_required = Math.max(0.8 * d_n_bolt * (Rp1 / cylinder_material[0]),
		0.8 * d_n_bolt); //11.4.3.3
	var F_bolt = F_pull / n_bolts;

	if(tube_end == "flanged" || tube_end == "trunnion") {
		F_bolt = Math.max(F_pull, Pa)/n_bolts;
	}

	var nominal_bolt_stress = Math.min(Rp1 / 3, Rm / 4);
	var A_s_bolt_required = F_bolt / nominal_bolt_stress;
	var E_n_min_bolt = d_n_bolt - 0.6495 * L_p_bolt;
	// TODO: Check if yield at design temp. should be added
	var f_thread_bolt = Math.min(cylinder_material[0] / 1.5, 
		cylinder_material[2] / 2.4, Rp1 / 1.5, Rm / 2.4);
	var f_sh_bolt = 0.8 * f_thread_bolt;
	var tau_shear_bolt = (Rp1 * A_s_bolt * 0.7 + 0.25 * F_bolt) 
		/ (0.5 * l_thread_bolt * E_n_min_bolt * Math.PI);

	set(part_name + "_" + "l_thread_bolt_required", l_thread_bolt_required);
	set(part_name + "_" + "A_s_bolt_required", A_s_bolt_required);
	set(part_name + "_" + "tau_shear_bolt", tau_shear_bolt);
	set(part_name + "_" + "f_sh_bolt", f_sh_bolt);

	l_thread_bolt_required= roundToDecimal(l_thread_bolt_required, 3);
	tau_shear_bolt = roundToDecimal(tau_shear_bolt, 3);
	f_sh_bolt= roundToDecimal(f_sh_bolt, 3);
	A_s_bolt_required= roundToDecimal(A_s_bolt_required, 3);

	if (l_thread_bolt_required > l_thread_bolt && tau_shear_bolt > f_sh_bolt) {
		error([], "The required engagement length of the " + part_name 
		+ " bolts are " + l_thread_bolt_required + " mm > " + l_thread_bolt 
		+ " mm and shear stress in the bolt is higher than the allowed \
		stress: " + tau_shear_bolt + " MPa > " + f_sh_bolt + " MPa. At least \
		one of these requirements has to be fulfilled. Refer to EN 13445-3 \
		11.4.3.3 and EN 14359 5.5.6.3, respectfully.");
	} else if (l_thread_bolt_required < l_thread_bolt && tau_shear_bolt < f_sh_bolt) {
		warn([], "The required engagement length of the " + part_name 
		+ " bolts is lower than or equal the actual length: " + l_thread_bolt_required 
		+ " mm =< " + l_thread_bolt + " mm" + " and shear stress in the bolt is \
		lower than or equal the critical stress: " + tau_shear_bolt + " MPa =< " 
		+ f_sh_bolt + " MPa. Refer to EN 13445-3 11.4.3.3 and EN 14359 5.5.6.3, respectfully.");
	} else if (l_thread_bolt_required > l_thread_bolt && tau_shear_bolt < f_sh_bolt) {
		warn([], "The required engagement length of the " + part_name 
		+ "  bolts is lower than or equal the actual length: " + l_thread_bolt_required 
		+ " mm =< " + l_thread_bolt + " mm. Refer to EN 13445-3 11.4.3.3 and EN 14359 5.5.6.3, \
		respectfully.");
	} else {
		warn([], "The shear stress in the bolt is lower than or equal the critical stress: " 
		+ tau_shear_bolt + " MPa =< " + f_sh_bolt + " MPa. Refer to EN 13445-3 11.4.3.3 \
		and EN 14359 5.5.6.3, respectfully.");
	}

	if (A_s_bolt_required > A_s_bolt) {
		error([], "The required stress area of the bolts is higher than the actual area: " 
		+ A_s_bolt_required + " mm^2 > " + A_s_bolt + " mm^2. Please refer to EN 13445-3 \
		11.4.3.3.");
	} else {
		warn([], "The required stress area of the bolts is lower or equal to the actual \
		area: " + A_s_bolt_required + " mm^2 =< " + A_s_bolt + " mm^2. Please refer to \
		EN 13445-3 11.4.3.3.");
	}
} // End of calc_bolts function

function bucklingStatic(cylinder) {
	var pi = Math.PI;
	var tube = cylinder.tube;
	var rod = cylinder.rod;
	var piston = cylinder.piston;
	var safetyFactor = 0;
	if (cylinder.bucklingMethod == "simple_static"){
		safetyFactor = safetyFactorSimple(cylinder.force, cylinder.maxLength/cylinder.euler, 
			tube.length/tube.euler, tube.inertia, rod.maxLength/rod.euler, rod.inertia, undefined);
	} else if (cylinder.bucklingMethod == "accurate_static") {
		safetyFactor = safetyFactorAccurate(cylinder.force, cylinder.maxLength/cylinder.euler, 
			cylinder.mass, tube.length/tube.euler, tube.inertia, rod.maxLength/rod.euler,
			rod.area, rod.inertia, rod.yield, rod.outerD, rod.end.innerD, rod.end.friction, 
			piston.guidingL/cylinder.euler, undefined, undefined, pi*rod.maxLength/cylinder.maxLength);
	} else if (cylinder.bucklingMethod == "en_static") {
		safetyFactor = safetyFactorEN(cylinder.force, cylinder.maxLength/cylinder.euler, 
			tube.length/tube.euler, tube.area, tube.inertia, tube.yield, rod.maxLength/rod.euler,
			rod.area, rod.inertia, rod.yield, undefined, undefined);
	}
	return safetyFactor;
}

function bucklingForceCurve(cylinder) {
	var pi = Math.PI;
	var tube = cylinder.tube;
	var rod = cylinder.rod;
	var piston = cylinder.piston;
	var forceCurve = cylinder.forceCurve;
	var stroke = undefined;
	var force = undefined;
	var cylinderLength = undefined;
	var rodLength = undefined;
	var guidingLength = undefined;
	var data = [];
	var safetyFactor = undefined;
	if (forceCurve.length < 2) {
		return {safetyFactor: undefined, stroke: undefined, force: undefined};
	}
	for (var i = 0; i < forceCurve.length; i++) {
		stroke = forceCurve[i].stroke;	
		force = forceCurve[i].f * Math.pow(10, 3); // Force in N
		cylinderLength = cylinder.minLength + stroke;
		rodLength = rod.minLength + stroke;
		guidingLength = rod.maxLength - rodLength + piston.guidingL;
		if (cylinder.bucklingMethod == "simple_force") {
			safetyFactor = safetyFactorSimple(force, cylinderLength/cylinder.euler, 
				tube.length/tube.euler, tube.inertia, rodLength/rod.euler, 
				rod.inertia, undefined);
		} else if (cylinder.bucklingMethod == "accurate_force") {
			safetyFactor = safetyFactorAccurate(force, cylinderLength/cylinder.euler, 
				cylinder.mass, tube.length/tube.euler, tube.inertia, rodLength/rod.euler,
				rod.area, rod.inertia, rod.yield, rod.outerD, rod.end.innerD, rod.end.friction, 
				guidingLength/cylinder.euler, undefined, undefined, pi*rodLength/cylinderLength);
		} else if (cylinder.bucklingMethod == "en_force") {
			safetyFactor = safetyFactorEN(force, cylinderLength/cylinder.euler, 
				tube.length/tube.euler, tube.area, tube.inertia, tube.yield, rodLength/rod.euler,
				rod.area, rod.inertia, rod.yield, undefined, undefined);
		}
		data.push({safetyFactor: safetyFactor, stroke: stroke, force: force});
	}
	data.sort((a,b) => (a.safetyFactor >= b.safetyFactor) ? 1 
		: (a.safetyFactor < b.safetyFactor) ? -1 : 0);
	return data[0];
}

function bucklingPressureCurve(cylinder) {
	var pi = Math.PI;
	var tube = cylinder.tube;
	var rod = cylinder.rod;
	var piston = cylinder.piston;
	var pressureCurve = cylinder.pressureCurve;
	var stroke = undefined;
	var force = undefined;
	var cylinderLength = undefined;
	var rodLength = undefined;
	var guidingLength = undefined;
	var data = [];
	var safetyFactor = undefined;
	if (pressureCurve.length < 2) {
		return {safetyFactor: undefined, stroke: undefined, force: undefined};
	}
	for (var i = 0; i < pressureCurve.length; i++) {
		stroke = pressureCurve[i].stroke;
		force = pressureCurve[i].p * piston.area * Math.pow(10, -1); // force in N
		rodLength = rod.minLength + stroke;
		cylinderLength = tube.length + rodLength;
		guidingLength = rod.maxLength - rodLength + piston.guidingL;
		if (cylinder.bucklingMethod == "simple_pressure") {
			safetyFactor = safetyFactorSimple(force, cylinderLength/cylinder.euler, 
				tube.length/tube.euler, tube.inertia, rodLength/rod.euler, 
				rod.inertia, undefined);
		} else if (cylinder.bucklingMethod == "accurate_pressure") {
			safetyFactor = safetyFactorAccurate(force, cylinderLength/cylinder.euler, 
				cylinder.mass, tube.length/tube.euler, tube.inertia, rodLength/rod.euler,
				rod.area, rod.inertia, rod.yield, rod.outerD, rod.end.innerD, rod.end.friction, guidingLength/cylinder.euler,
				undefined, undefined, pi*rodLength/cylinderLength);
		} else if (cylinder.bucklingMethod == "en_pressure") {
			safetyFactor = safetyFactorEN(force, cylinderLength/cylinder.euler, 
				tube.length/tube.euler, tube.area, tube.inertia, tube.yield, rodLength/rod.euler,
				rod.area, rod.inertia, rod.yield, undefined, undefined);
		}
		data.push({safetyFactor: safetyFactor, stroke: stroke, force: force});
	}
	data.sort((a,b) => (a.safetyFactor >= b.safetyFactor) ? 1 
		: (a.safetyFactor < b.safetyFactor) ? -1 : 0);
	return data[0];
}

function safetyFactorSimple(actualLoad, cylinderLength, tubeLength, tubeInertia, rodLength,
	rodInertia, E) {
	/* Calculates the buckling safety factor according to the method in [3.2.2]. */
	if (typeof E == "undefined") {
		E = 206000;
	}

	var bucklingLoad = calcBucklingLoadSimple(tubeLength, tubeInertia, rodLength, rodInertia,
		cylinderLength, E);
	var safetyFactor = bucklingLoad/actualLoad;
  console.log("------------- In safetyFactorSimple -------------")
  console.log("- actual load: ", actualLoad)
  console.log("- buckling load: ", bucklingLoad)
  console.log("- safetyFactor: ", safetyFactor)
  console.log("- L: ", cylinderLength)
  console.log("- L1: ", tubeLength)
  console.log("- L2: ", rodLength)
  console.log("- I1: ", tubeInertia)
  console.log("- I2: ", rodInertia)
  console.log("- E: ", E)
  return roundToDecimal(safetyFactor, 3);
}

function calcBucklingLoadSimple(L1, I1, L2, I2, L, E) {
	/* Calculates the critical buckling load of a cylinder based the method 
	 * in DNVGL-ST-0194 [3.2.2].
	 * arg L1: float, tube length [mm]
	 * arg I1: float, tube inertia moment [mm^4]
	 * arg L2: float, rod length [mm]
	 * arg I2: float, rod inertia moment [mm^4]
	 * arg L: float, cylinder length between bearings [mm]
	 * arg E: float, Young's modulus */
	var pi = Math.PI;
	var Z = L1/I1 + L2/I2 + (1/I2 - 1/I1)*L/(2*pi)*Math.sin(2*pi*L1/L);
	var Fe = E*pi*pi/(L*Z);
	return Fe;
}

function safetyFactorAccurate(actualLoad, cylinderLength, cylinderMass, tubeLength, tubeInertia, 
	rodLength, rodArea, rodInertia, rodYield, rodDiameter, rodEndEyeDiameter, frictionCoefficient, pistonGuideLength, 
	E, delta, alpha) {
	/* Calculates the buckling safety factor according to the method in [A.4.3]. */
	var rodEndEyeRadius = null;
	if (typeof E == "undefined") {
		E = 206000;
	}
	if (typeof delta == "undefined") {
		delta = 0.17;
	}
	if (typeof rodEndEyeDiameter == "undefined") {
		rodEndEyeRadius = 0;
	} else {
		rodEndEyeRadius = rodEndEyeDiameter/2;
	}

	var bucklingLoad = calcBucklingLoadAccurate(cylinderLength, tubeLength, tubeInertia,
		rodLength, rodInertia, pistonGuideLength, rodArea, rodYield, rodDiameter,
		rodEndEyeRadius, frictionCoefficient, cylinderMass, E, delta, alpha);
	var safetyFactor = bucklingLoad/actualLoad;
		console.log("In safetyFactorAccurate: ")
		console.log("- bucklingLoad: ", bucklingLoad)
		console.log("- actualLoad: ", actualLoad)
		console.log("- safetyFactor: ", safetyFactor)
	return roundToDecimal(safetyFactor, 3);
}

function calcBucklingLoadAccurate(L, L1, I1, L2, I2, L5, A, R, d, r, my, m, E, delta, alpha) {
	/* Calculates the buckling load according to [A.4.3].
	 * arg L: float, cylinder length [mm]
	 * arg L1: float, tube length [mm]
	 * arg I1: float, tube inertia [mm^4]
	 * arg L2: float, rod length [mm]
	 * arg I2: float, rod inertia [mm^4] 
	 * arg L5: float, piston guide length [mm]
	 * arg A: float, rod area [mm^2]
	 * arg R: float, rod yield stress [MPA]
	 * arg d: float, rod outer diameter [mm]
	 * arg r: float, rod end eye inner diameter [mm]
	 * arg my: float, rod eye friction coefficient [-]
	 * arg m: float, cylinder mass [kg]
	 * arg E: float, E-modulus [MPA]
	 * arg delta: float, clearance [-]
	 * return: float*/
	var pi = Math.PI;
	var g = 9.81;

	var AA = L1/(2*I1) + L2/(2*I2) + L/(4*pi) * (1/I1 - 1/I2) * Math.sin(2*alpha);
	var BB = 4*L/(3*pi) * (1/I2 - 1/I1) * Math.pow(Math.sin(alpha), 3);
	var CC = L1/(2*I1) + L2/(2*I2) + L/(8*pi)*(1/I1 - 1/I2)*Math.sin(4*alpha);
	var DD = Math.pow(pi, 2)*E/(2*L);

	var a = 4*AA*CC-BB*BB;
	var b = -4*DD*(4*AA+CC);
	var c = 16*DD*DD;

	var Pe = (-b-Math.sqrt(b*b-4*a*c))/(2*a);
	var f = L1*delta/(2*L*L5) * Math.sqrt(pi*pi*E*I2/Pe)

	var FF = d/(2*I2)*(my*r + f);
	var GG = m*g*L*d/(16*I2);
	var HH = Math.sqrt(1 + 2*A*FF - 2*R*A/Pe + A*A*FF*FF + 2*A*A*FF*R/Pe
		+ Math.pow(R*A/Pe, 2) + 4*A*GG/Pe);
	var P = R*A/2 + (Pe/2)*(1 + A*FF - HH);
		console.log("------------- In calcBucklingLoadAccurate -------------")
		console.log("- L: ", L)
		console.log("- L1: ", L1)
		console.log("- L2: ", L2)
		console.log("- L5: ", L5)
		console.log("- A: ", A)
		console.log("- r: ", r)
		console.log("- d: ", d)
		console.log("- Reh: ", R)
		console.log("- my: ", my)
		console.log("- m: ", m)
		console.log("- E: ", E)
		console.log("- alpha: ", alpha)
		console.log("- delta: ", delta)
		console.log("- AA: ", AA)
		console.log("- BB: ", BB)
		console.log("- CC: ", CC)
		console.log("- DD: ", DD)
		console.log("- a: ", a)
		console.log("- b: ", b)
		console.log("- c: ", c)
		console.log("- Pe: ", Pe)
		console.log("- f: ", f)
		console.log("- FF: ", FF)
		console.log("- GG: ", GG)
		console.log("- HH: ", HH)
		console.log("- P: ", P)
	return P;
}

function safetyFactorEN(actualLoad, cylinderLength, tubeLength, tubeArea, tubeInertia, tubeYield, 
	rodLength, rodArea, rodInertia, rodYield, E, alpha) {
	/* Calculates the buckling safety factor of a hydraulic cylinder based on the 
	 * buckling curve from EN 1993-1-1 as referred to in DNVGL-ST-0194 [A.5].
	 * arg actualLoad: float, actual maximum pushing force [N]
	 * arg cylinderLength: float, cylinder length between bearings [mm]
	 * arg tubeArea: float, tube cross section area [mm^2]
	 * arg tubeYield: float, tube yield stress [MPa]
	 * arg rodArea: float, rod cross section area [mm^2]
	 * arg rodYield: float, rod yield stress [MPA]
	 * arg alpha: imperfection factor [-]
	 * return: float */
	if (typeof E == "undefined") {
		E = 206000;
	}
	if (typeof alpha == "undefined") {
		alpha = 0.49;
	}
  console.log("------------- In safetyFactorEN -------------")
	var criticalLoad = calcBucklingLoadSimple(tubeLength, tubeInertia, rodLength, rodInertia, 
		cylinderLength, E);
  console.log("- E: ", E)
  console.log("- alpha: ", alpha)
  console.log("- tubeYield: ", tubeYield)
  console.log("- rodYield: ", rodYield)  
  console.log("- actual load: ", actualLoad)
	console.log("- critical load: ", criticalLoad)
  console.log("- For tube: ")
	var tubeCapacity = calcBucklingCapacity(criticalLoad, tubeArea, tubeYield, alpha);
	var tubeLoad = tubeCapacity * tubeArea * tubeYield;
  console.log("- tubeLoad: ", tubeLoad)
  console.log("- For rod: ")
	var rodCapacity = calcBucklingCapacity(criticalLoad, rodArea, rodYield, alpha);
	var rodLoad = rodCapacity * rodArea * rodYield;
  console.log("- rodLoad: ", rodLoad)
	var bucklingLoad = Math.min(tubeLoad, rodLoad);
	var safetyFactor = bucklingLoad/actualLoad;		
  console.log("- buckling load: ", bucklingLoad)
  console.log("- safety factor: ", safetyFactor)
	return roundToDecimal(safetyFactor, 3);
}

function calcBucklingCapacity(Fe, A, Reh, alpha) {
	/* Calculates the buckling capacity of a circular part based on the method in
	 * DNVGL-ST-0194 [A.5.3].
	 * arg Fe: float, buckling load [N]
	 * arg A: float, cross section area [mm^2]
	 * arg Reh: float, yield strength [MPA]
	 * arg alpha: float, imperfection factor (0.49 for buckling curve c) [-] */
	var lambda = Math.sqrt(A*Reh/Fe);
	var phi = 0.5 * (1 + alpha*(lambda - 0.2) + Math.pow(lambda,2));
	var chi = 1 / (phi + Math.sqrt(Math.pow(phi,2) - Math.pow(lambda,2)));
  console.log("   - lambda: ", lambda)
  console.log("   - phi: ", phi)
  console.log("   - chi: ", chi)
	return chi;
}



////////////////////////////////////////////////////////////////
//////// GENERAL SHORTCUTS FOR CREATING PARTS OF THE GUI ///////
////////////////////////////////////////////////////////////////

var uiMaterial = {
  type: "dropdown",
  options: function options(get, set) {
    if (get("material_define") == "standard") {
      var items = [{ text: "S355J2+N, EN10025-3",
        value: {
          mat: "S355J2+N, EN10025-3"
        }
      }, {
        text: "E355+SR, EN10305-1",
        value: {
          mat: "E355+SR, EN10305-1"
        } }, {
        text: "1.4462, EN 10088-3",
        value: {
          mat: "1.4462, EN 10088-3"
        }
      }, {
        text: "1.4460, EN 10088-3",
        value: {
          mat: "1.4460, EN 10088-3"
        }
      }, {
        text: "1.4418 +QT760 EN 10088-3",
        value: {
          mat: "1.4418 +QT760 EN 10088-3"
        }
      }, {
        text: "1.4418 +QT900 EN 10088-3",
        value: {
          mat: "1.4418 +QT900 EN 10088-3"
        }
      }, {
        text: "1.4404 EN 10088-3",
        value: {
          mat: "1.4404 EN 10088-3"
        }
      }, {
        text: "1.4403 EN 10088-3",
        value: {
          mat: "1.4403 EN 10088-3"
        }
      }];
    } else {
      var items = [{
        text: "Carbon steel",
        value: {
          mat: "cs"
        }
      }, {
        text: "Austenetic stainless steel",
        value: {
          mat: "ss"
        }
      }];
    }
    return items;
  }
} // End of uiMaterial dropdown menu type

var uiBoltMaterial = {
  type: "dropdown",
  options: function options(get, set) {
    var items = [{ text: "8.8", value: { mat: "8.8", f_y: 640, TS: 800 } }, 
	    { text: "10.9", value: { mat: "10.9", f_y: 900, TS: 1040 } }];
    
    return items;
  }
} // End of uiBoltMaterial dropdown menu type

var uiCalcMethod = {
  type: "dropdown",
  options: function options(get) {
    var items = [{ text: "Method [3.3.2]", value: "simple_static" }, 
	    { text: "Method [3.3.2], stroke-force curve", value: "simple_force" }, 
	    { text: "Method [3.3.3], stroke-pressure curve", value: "simple_pressure" },
	    { text: "Method [A.4]", value: "accurate_static" },
	    { text: "Method [A.4], stroke-force curve", value: "accurate_force" }, 
	    { text: "Method [A.4], stroke-pressure curve", value: "accurate_pressure" }, 
	    { text: "Method [A.5]", value: "en_static" }, 
	    { text: "Method [A.5], stroke-force curve", value: "en_force" }, 
	    { text: "Method [A.5], stroke-pressure curve", value: "en_pressure" }];
    
    return items;
  }
} // End of uiCalcMethod dropdown menu type

var uiListMaterial = {
  type: "dropdown",
  options: function options(get) {
    var mat_name = get("mat_name"),
        Rp1 = get("Rp1"),
	Rp1t = get("Rp1t"),
        Reh = get("Reh"),
	Ret = get("Ret"),
        Rm = get("Rm"),
        mat = get("material_type"),
        def = get("material_define");

    if (!Array.isArray(mat_name)) {
      mat_name = [mat_name];
      Rp1 = [Rp1];
      Rp1t = [Rp1t];	
      Reh = [Reh];
      Ret = [Ret];
      Rm = [Rm];
      mat = [mat];
      def = [def];
    }
    var items = [];
    for (var i = 0; i < mat_name.length; i++) {
      var mat_type;

      if (mat[i] == undefined)
      {
        continue;
      }

      if (def[i] == "standard") {
        mat_name[i] = mat[i].mat;
        if (mat_name[i] == MatList[0]) {
          var f_y = f_y_S_355;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[1]) {
          var f_y = f_y_E_355;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[2]) {
          var f_y = f_y_1_4462;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[3]) {
          var f_y = f_y_1_4460;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[4]) {
          var f_y = f_y_1_4418QT760;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[5]) {
          var f_y = f_y_1_4418QT900;
          mat_type = "cs";
        } else if (mat_name[i] == MatList[6]) {
          var f_y = f_y_1_4404;
          mat_type = "ss";
        } else {
          mat_type = "ss";
          var f_y = f_y_1_4403;
        }
      } else {
        if (mat[i].mat == "ss") {
          var f_y = [[3, 400, Rp1[i], Rp1t[i], Rm[i]]];
          mat_type = mat[i].mat;
        } else {
          var f_y = [[3, 400, Reh[i], Ret[i], Rm[i]]];
          mat_type = mat[i].mat;
        }
      }

      items.push({
        text: mat_name[i],
        value: {
          f_y: f_y,
          Rp1: Rp1[i],
	  Rp1t: Rp1t[i],
          Reh: Reh[i],
	  Ret: Ret[i],
          Rm: Rm[i],
          mat: mat_type,
          name: mat_name[i],
	  def: def[i]
        }
      });
    }
    return items;
  }
} // End of uiListMaterial dropdown menu type

function itText(name, desc, v, help, show, img) {
  return {
    name: name,
    type: "input",
    ui: "text",
    desc: desc,
    var: v,
    show: show,
    help: help == undefined ? undefined : { text: help, img: img }
  };
} // End of itText function

function itNumPos(name, desc, v, help, show, img) {
  return {
    name: name,
    type: "input",
    ui: {
      "type": "number",
      "min": 0
    },
    desc: desc,
    var: v,
    show: show,
    help: help == undefined ? undefined : { text: help, img: img }
  };
} // End of itNumPos function

function itNum(name, desc, v, help, show, img) {
  return {
    name: name,
    type: "input",
    ui: {
      "type": "number"
    },
    desc: desc,
    var: v,
    show: show,
    help: help == undefined ? undefined : { text: help, img: img }
  };
} // ENd of itNum function

////////////////////////////////////////////////////////////////
///////////////////// APPLICATION FUNCTION /////////////////////
////////////////////////////////////////////////////////////////

define(function () {
  
  //////////////////////////////////////////////////////////////////////////////
  ////// DEFNINITION OF THE GRAPHICAL USER INTERFACE THAT SHOULD BE SHOWN //////
  //////////////////////////////////////////////////////////////////////////////

  return {
    title: "Hydraulic Cylinders, CG-0194",
    edition: "December, 2015",
    desc: "Calculation, submital and approval tool for hydraulic cylinders",
    intro: '\nWelcome to the PASS - Pre-Approval Screening Sevice.\n\nYou may use this calculation tool to check your design against DNV GL requirements and receive an early warnings of possible non-conformities before submission to design review.\n',
    app:
    {
      files: {
        icon: 'icon',
        menuIcon: 'menuIcon'
      },
      desc: "A simple calculation tool to allow you to receive an instant yes/no check of your hydraulic cylinder design as part of the drawing submittal process",
      interface: 
      {
        // Interface menu blocks
        blocks:
        [
          // Project information menu block
          {
            title: "Project information",
            name: "info",
            items: 
            [
              // Type/Model input block
              itText("model", "Type/model", undefined, 
                "Please enter the hydraulic cylinder identifying model or tag number"),
              {
                name: "forship",
                type: "input",
                desc: "For DNV GL classed vessel?",
                ui: {
                  type: "dropdown",
                  options: [{ text: "Yes", value: "yes" }, { text: "No", value: "no" }]
                }
              },
              // Vessel ID input block
              itText("shipid", "DNV GL vessel id", undefined, 
                "The vessel identification number is usually in the form of D followed by \
                five numbers for newbuild projects.",
                function (get) {return get("forship") == "yes";}),
              // Manufacturer input block
              itText("manufacturer", "Manufacturer", undefined, 
                "Please enter your business name or DNV GL customer number"),
              // Hydraulic cylinder application input block
              {
                name: "pv_func",
                type: "input",
                desc: "Hydraulic cylinder application",
                help: {
                  text: "The identifying function corresponds to the DNV GL DocReq and helps us to \
                    ensure we have fulfilled the documentation requirements for the project."
                },
                ui: {
                  type: "dropdown",
                  options: [ {
                    text: "Hydraulic cylinders for general machinery use",
                    value: "machinery"
                  },{
                    text: "Hydraulic cylinders for steering gear/water jet steering",
                    value: "waterjet"
                  }, {
                    text: "Hydraulic cylinders for jacking systems for self elevating units",
                    value: "self_elevating"
                  }, {
                    text: "Hydraulic cylinders for offshore and platform cranes",
                    value: "offshore"
                  }, {
                    text: "Hydraulic cylinders for water tight doors (internal locking mechanism)",
                    value: "water_tight_doors"
                  }, {
                    text: "Hydraulic cylinders for offshore gangways",
                    value: "gangway"
                  }]
                }
              },
              // GA drawing number input block
              itText("dwg", "GA drawing number", undefined, 
                "Please give a reference to the main drawing that will be submitted for approval"),
              // GA revision input block
              itText("dwgrev", "Revision of GA drawing", undefined, 
                "Please give a reference to the main drawing that will be submitted for approval"),
              // Manufacturer reference input block
              itText("project", "Manufacturer reference", undefined, "Your project reference")
            ]
            // End of project information menu items block
          },
          // End of project information menu block

          // Design data menu block
          {
            title: "Design data",
            name: "basic",
            desc: "The hydraulic cylinder design parameters",
            items: 
            [
              // Design pressure push input block
              itNumPos("DPpush", "Design pressure, push (bar)", "P_{d.push}"), 
              // Design pressure pull input block
              itNumPos("DPpull", "Design pressure, pull (bar)", "P_{d.pull}"), 
              // Test pressure push input block
              itNumPos("p_test_push", "Test pressure, push (bar)", "P_{d.push}"), 
              // Test pressure pull input block
              itNumPos("p_test_pull", "Test pressure, pull (bar)", "P_{d.pull}"), 
              // Minimum design temperature input block
              itNum("DTmin", "Minimum design temperature (degC)", "T_{d.min}", 
                "Charpy v-notch impact testing is required at the minimum design temperature. \
                Average 27 joules transverse is required. Impact testing is not required for \
                steering cylinders."), 
              // Maximum design temperature input block
              itNum("DT", "Maximum design temperature (degC)", "T_d"),
              // Corrosion allowance output block
              {
                name: "c",
                type: "output",
                desc: "Corrosion allowance (mm)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  return 0.3
                },
                help: {
                  text: "For carbon and low-alloy steels, the corrosion allowance c shall be in \
                    general 0.3 mm. For austenitic steels or titanium no corrosion allowance will \
                    in general be required. Where adverse corrosion or erosion conditions exist a \
                    greater allowance may be required. For pressurised parts with more than 30 mm \
                    material thickness made from carbon or carbon manganese steels the corrosion \
                    allowance can be exempted."
                }
              },
              // Medium type dropdown block
              {
                name: "medium_type",
                type: "input",
                desc: "Medium type",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Hydraulic fluid, oil based",
                    value: "Hydraulic fluid, oil based"
                  }, {
                    text: "Hydraulic fluid, non-flammable water based",
                    value: "Hydraulic fluid, non-flammable water based"
                  }, {
                    text: "Pneumatic, air",
                    value: "air"
                  }]
                }
              },
              // Working pressure input block
              itNumPos("wp", "Working pressure (bar)", "P_w", undefined, 
                function (get) {return get("pv_func") == "waterjet";})
            ]
            // End of design data menu items block
          },
          // End of design data menu block

          // Material data menu block
          {
		title: "Material data",
		name: "mblock",
		intro: "test",
		multiple: true,
		desc: "Add your materials for the project",
		items: 
            [ 
              // Material standard dropdown block
		{ 
			name: "material_define",
			type: "input",
			desc: "Material custom/standard",
			help: {
			  text: "Materials for cylinder tube, piston rod, end covers and end eyes shall \
			    be delivered with 3.1 material certificates (according to EN 10204 or equivalent). \
			    Approval of manufacturer is not required. For cylinders intended for steering gear \
			    or water jet steering applications, 3.2 material certification is required for the \
			    cylinder tube and piston rod. The end covers, end eye and piston shall be delivered \
			    with 3.1 certificates."
			},
			ui: {
			  type: "dropdown",
			  options: [{
			    text: "Standard materials (at room temperature)",
			    value: "standard"
			  }, {
			    text: "Custom defined materials",
			    value: "custom"
			  }]
			}
		},
              // Material type dropdown block
              {
                name: "material_type",
                type: "input",
                desc: "Material type",
                help: {
                  text: "Materials for cylinder tube, piston rod, end covers and end eyes shall be \
                    delivered with 3.1 material certificates (according to EN 10204 or equivalent). \
                    Approval of manufacturer is not required. For cylinders intended for steering \
                    gear or water jet steering applications, 3.2 material certification is required \
                    for the cylinder tube and piston rod. The end covers, end eye and piston shall \
                    be delivered with 3.1 certificates."
                },
                ui: uiMaterial
              },
              itText("mat_name", "Material grade and standard", undefined, undefined, 
                function (get) {return get("material_define") == "custom";}),
              itNumPos("Rp1", "1% proof stress at room temperature (MPa)", "R_{p1.0}", undefined, 
                function (get) {return get("material_define") == "custom" && get("material_type").mat == "ss";}),
              itNumPos("Reh", "Yield stress at room temperature (MPa)", "R_{eH}", undefined, 
                function (get) {return get("material_define") == "custom" && get("material_type").mat == "cs";}), 
	      itNumPos("Rp1t", "1% proof stress at design temperature (MPa)", "R_{p1.0t}", undefined, 
                function (get) {return get("material_define") == "custom" && get("material_type").mat == "ss";}),
              itNumPos("Ret", "Yield stress at design temperature (MPa)", "R_{et}", undefined, 
                function (get) {return get("material_define") == "custom" && get("material_type").mat == "cs";}),
              itNumPos("Rm", "Tensile strength at room temperature (MPa)", "R_{m}", undefined, 
                function (get) {return get("material_define") == "custom";})
            ]
            // End of material data menu items block
          }, 
          // End of material data menu block

          // Dimensions menu block
          {
            title: "Dimensions",
            name: "stroke",
            desc: "Information about the cylinder's minimum length and maximum stroke",
            items: 
            [
              // End connection at tube end dropdown block
              {
                name: "tubeEnd",
                type: "input",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "End eye",
                    value: "endEye"
                  }, {
                    text: "Trunnion mounted",
                    value: "trunnion"
                  }, {
                    text: "Flanged (Located at cylinder tube)",
                    value: "flanged"
                  }, {
                    text: "Flanged (Located at end cover)",
                    value: "flanged_located_at_end_cover"
                  }]
                },
                desc: "End connection at the tube end",
                help: {
                  text: "Attachment at the tube end",
                  img: "ends.png"
                }
              },
              // End connection at rod end dropdown block
              {
                name: "rodEnd",
                type: "input",
                ui: {
                  type: "dropdown",
                  options: 
                  [
                    {
                      text: "End eye",
                      value: "endEye"
                    },
                    {
                      text: "Threaded (fixed)",
                      value: "threaded_fixed"
                    },
                    {
                      text: "Threaded (pinned)",
                      value: "threaded_pinned"
                    }
                  ]
                },
                desc: "End connection at the rod end",
                help: {
                  text: "Attachment at the rod end"
                }
              },
              // Rod always inside menu for "Trunnion mounted" cylinders
              {
                name: "rodInside",
                type: "input",
                show: function(get){return (get("tubeEnd") == "trunnion");},
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "No",
                    value: "no"
                  }, {
                    text: "Yes",
                    value: "yes"
                  }]
                },
                desc: "Rod always inside cylinder"
              },
              // Cylinder figure block 
              { // Cylinder figure 1 (flanged, eye)
                name: "measuresfig1",
                type: "figure",
                src: flanged_endEye,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "flanged") && get("rodEnd") == "endEye";
                }
              }, { // Cylinder figure 2 (eye, eye)
                name: "measuresfig2",
                type: "figure",
                src: endEye_endEye,
                desc: "Cylinder length",
                show: function show(get) {
                  return get("tubeEnd") == "endEye" && get("rodEnd") == "endEye";
                }
              }, { // Cylinder figure 3 (trunnion, eye, outside)
                name: "measuresfig3",
                type: "figure",
                src: trunnion_endEye,
                desc: "Cylinder length",
                show: function show(get) {
                  return get("tubeEnd") == "trunnion" && get("rodEnd") == "endEye" && get("rodInside") == "no" ;
                }
              }, { // Cylinder figure 4 (trunnion, threaded, outside)
                name: "measuresfig4",
                type: "figure",
                src: trunnion_thread,
                desc: "Cylinder length",
                show: function show(get) {
                  return get("tubeEnd") == "trunnion" 
                    && (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned") 
                    && get("rodInside") == "no";
                }
              }, { // Cylinder figure 5 (eye, threaded)
                name: "measuresfig5",
                type: "figure",
                src: endEye_thread,
                desc: "Cylinder length",
                show: function show(get) {
                  return get("tubeEnd") == "endEye" 
                    && (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned");
                }
              }, { // Cylinder figure 6 (flanged, threaded)
                name: "measuresfig6",
                type: "figure",
                src: flanged_thread,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "flanged") 
                    && (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned");
                }
              }, { // Cylinder figure 7 (flanged end cover, eye)
                name: "measuresfig7",
                type: "figure",
                src: flanged_located_at_end_cover_endEye,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "flanged_located_at_end_cover") 
                    && get("rodEnd") == "endEye";
                }
              }, { // Cylinder figure 8 (flanged end cover, threaded)
                name: "measuresfig8",
                type: "figure",
                src: flanged_located_at_end_cover_thread,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "flanged_located_at_end_cover") 
                    && (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned");
                }
              }, { // Cylinder figure 9 (trunnion, eye, inside)
                name: "measuresfig9",
                type: "figure",
                src: flanged_rod_always_inside_cylinder,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "trunnion") 
                    && get("rodEnd") == "endEye" 
                    && get("rodInside") == "yes";
                }
              }, { // Cylinder figure 10 (trunnion, threaded, inside)
                name: "measuresfig10",
                type: "figure",
                src: flanged_rod_always_inside_cylinder_thread,
                desc: "Cylinder length",
                show: function show(get) {
                  return (get("tubeEnd") == "trunnion") 
                    && (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned") 
                    && get("rodInside") == "yes";
                }
              },
              // Length of cylinder part input block
              itNumPos("L1", "Length of cylinder part (mm)", "L_1",
                "Length of the cylinder part from the centre of its mounting (mm)", undefined, undefined),
              // Visible length of piston rod input block 
              itNumPos("L2max", "Visible length of piston rod, fully extracted position (mm)", "L_{2}", 
                "Visible length of piston rod in fully retracted position from centre of its mounting (mm)", 
                undefined, undefined),
              // Maximum stroke input block 
              itNumPos("Stroke", "Maximum stroke (mm)"),
              // Minimum visible piston rod length block
              {
                name: "L2mino",
                var: "L_{4}",
                type: "output",
                desc: "Minimum visible piston rod length (mm)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  return 1 * get("L2max") - 1 * get("Stroke");
                }
              }, 
              // Maximum length of cylinder block
              { 
                name: "Lo",
                var: "L",
                type: "output",
                desc: "Maximum length of cylinder",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  return 1 * get("L1") + 1 * get("L2max");
                }
              },
              // Buckling calculation method list block
              {
                name: "calcMethod",
                type: "input",
                desc: "Buckling calculation method",
                help: {
                  text: "The simple calculation method calculates the buckling safety \
                    factor at the maximum pressure and stroke. Stroke-pressure curve and \
                    stroke-force curve calculates the buckling safety factor incrementally \
                    based on different pressures/forces at different stroke lengths"
                },
                ui: uiCalcMethod
              },
              //
              itNumPos("L3", "Guiding length of piston in cylinder tube (mm)", "L_3", 
                "The Guiding length of piston in cylinder tube as given on the general arrangement drawing", function show(get) {
			var calcMethod = get("calcMethod");
                	return calcMethod == "accurate_static" || calcMethod == "accurate_force" || calcMethod == "accurate_pressure";
              }),
              //
              itNumPos("m_cyl", "The weigth of the cylinder, kg", "m", undefined, function show(get) {
		      	var calcMethod = get("calcMethod");
                	return calcMethod == "accurate_static" || calcMethod == "accurate_force" || calcMethod == "accurate_pressure";
              }),
              //
              itNumPos("my_end_eye", "Coefficient of friction for end eyes [Default: 0.19]", undefined, undefined, function show(get) {
			var calcMethod = get("calcMethod");
                	return calcMethod == "accurate_static" || calcMethod == "accurate_force" || calcMethod == "accurate_pressure";
              }), {
                name: "fCurve",
                type: "input",
                desc: "Enter the stroke points and corresponding force for evaluation",
                show: function show(get) {
			var calcMethod = get("calcMethod");
                  	return calcMethod == "simple_force" || calcMethod == "accurate_force" || calcMethod == "en_force";
                },
                ui: {
                  type: "table",
                  columns: [{ name: "stroke", title: "Stroke (mm)" }, { name: "f", title: "Pushing force (kN)" }],
                  rows: -1
                }
              }, {
                name: "pCurve",
                type: "input",
                desc: "Enter the stroke points and corresponding pressure for evaluation",
                show: function show(get) {
			var calcMethod = get("calcMethod");
                  return calcMethod == "simple_pressure" || calcMethod == "accurate_pressure" || calcMethod == "en_pressure";
                },
                ui: {
                  type: "table",
                  columns: [{ name: "stroke", title: "Stroke (mm)" }, { name: "p", title: "Pushing pressure (bar)" }],
                  rows: -1
                }
              },
		// Cylinder tube heading block 
		{
			name: "head1",
			type: "heading",
			desc: "Cylinder tube"
		}, 
		// Cylinder tube minimum outer diameter input block
		itNumPos("OD", "Cylinder tube minimum outer diameter (mm)", "D_o", "The outer diameter of the cylinder tube as given on the general arrangement drawing"), 
		// Cylinder tube maximum outer diameter input block
		itNumPos("DI", "Cylinder tube maximum inner diameter (mm)", "D_i", "The inner diameter of the cylinder tube as given on the general arrangement drawing"),
		// Cylinder tube material block
		{
			name: "tubeMaterial",
			type: "input",
			desc: "Select cylinder tube material",
			ui: uiListMaterial
		},
		// Weld joint efficiency block 
		{
			name: "v",
			var: "v",
			type: "input",
			ui: {
				type: "dropdown",
				options: [{ text: "1.0 for seamless tubes or 100% NDT (RT or UT)", value: 1 }, 
					{ text: "0.85 for spot or random NDT (RT or UT)", value: 0.85 }, 
					{ text: "0.7 for visual inspection", value: 0.7 }]
			},
			desc: "Weld joint efficiency",
			help: { // Help
				text: "Refer to DNV GL Ship Pt.4 Ch.7 Sec.7 [4.1] for the extent of corresponding NDT"
			}
		},
		// Cylinder tolerance and fabrication allowance input block 
		itNumPos("cm_shell", "Tolerance and fabrication allowance (mm)", "c_m", "The calculated thicknesses according to the formulae are the minimum required. Minus tolerance on wall thickness and a possible reduction in wall thickness due to forming shall be added to the calculated thicknesses."), 
		// Cylinder 1% proof stress at room temperature
		{
			name: "R_p_tube",
			var: "R_{p1}",
			type: "output",
			show: function show(get) {
				return get("tubeMaterial").mat == "ss";
			},
			desc: "1% proof stress at room temperature of cylinder (MPa)",
			ui: {
				type: "text"
			},
			value: function value(get) {
				var thickness = (get("OD") - get("DI")) / 2;
				var materialType = get("tubeMaterial").name;
				var material = get("tubeMaterial");
				if (materialType == "E355+SR, EN10305-1") {
					var thickness = get("OD");
				}
				var data = getYieldAndStrength(thickness, material);
				return data[0];
			}
		}, 
		// Cylinder 1% proof stress at design temperature
		{
			name: "R_p1t_tube",
			var: "R_{p1t}",
			type: "output",
			show: function show(get) {
				var tubeMaterial = get("tubeMaterial");
				return tubeMaterial.mat == "ss" && tubeMaterial.def == "custom";
			},
			desc: "1% proof stress at design temperature of cylinder (MPa)",
			ui: {
				type: "text"
			},
			value: function value(get) {
				var thickness = (get("OD") - get("DI")) / 2;
				var materialType = get("tubeMaterial").name;
				var material = get("tubeMaterial");
				if (materialType == "E355+SR, EN10305-1") {
					var thickness = get("OD");
				}
				var data = getYieldAndStrength(thickness, material);
				return data[1];
			}
		},   
		// Cylinder yield stress at room temperature
		{
			name: "f_y_tube",
			var: "R_{eh}",
			type: "output",
			show: function show(get) {
				return get("tubeMaterial").mat == "cs";
			},
			desc: "Yield strength of the cylinder at room temperature (MPa)",
			ui: {
				type: "text"
			},
			value: function value(get) {
				var thickness = (get("OD") - get("DI")) / 2;
				var materialType = get("tubeMaterial").name;
				var material = get("tubeMaterial");
				if (materialType == "E355+SR, EN10305-1") {
					var thickness = get("OD");
				}
				var data = getYieldAndStrength(thickness, material);
				return data[0];
			}
		},
		    // Cylinder yield stress at design temperature
		{
			name: "f_yt_tube",
			var: "R_{et}",
			type: "output",
			show: function show(get) {
				var tubeMaterial = get("tubeMaterial");
				return tubeMaterial.mat == "cs" && tubeMaterial.def == "custom";
			},
			desc: "Yield strength of the cylinder at design temperature (MPa)",
			ui: {
				type: "text"
			},
			value: function value(get) {
				var thickness = (get("OD") - get("DI")) / 2;
				var materialType = get("tubeMaterial").name;
				var material = get("tubeMaterial");
				if (materialType == "E355+SR, EN10305-1") {
					var thickness = get("OD");
				}
				var data = getYieldAndStrength(thickness, material);
				return data[1];
			}
		},
              // Cylinder tensile strength output block
              {
                name: "TS_tube",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the cylinder at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD") - get("DI")) / 2;
                  var materialType = get("tubeMaterial").name;
                  var material = get("tubeMaterial");
                  if (materialType == "E355+SR, EN10305-1") {
                    var thickness = get("OD");
                  }
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // Rod heading block
              {
                name: "head2",
                type: "heading",
                desc: "Rod"
              }, 
              // Rod material list block
              {
                name: "rodMaterial",
                type: "input",
                desc: "Select rod material",
                ui: uiListMaterial
              },
              // Piston rod outer diameter input block
              itNumPos("OD_rod", "Piston rod outer diameter (mm)", "d_o"), 
              // Piston rod inner diameter input block
              itNumPos("DI_rod", "Piston rod inner diameter (mm)", "d_i", "Consider 0 mm for a solid rod"), 
              // Rod 1% proof stress at room temperature
              {
                name: "R_p_rod",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("rodMaterial").mat == "ss";
                },
                desc: "1% proof stress of the rod at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD_rod") - get("DI_rod")) / 2;
                  var materialType = get("rodMaterial").mat;
                  var material = get("rodMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
              // Rod 1% proof stress at design temperature
              {
                name: "R_p1t_rod",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var rodMaterial = get("rodMaterial");
                  	return rodMaterial.mat == "ss" && rodMaterial.def == "custom";
                },
                desc: "1% proof stress of the rod at design temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD_rod") - get("DI_rod")) / 2;
                  var materialType = get("rodMaterial").mat;
                  var material = get("rodMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[1];
                }
              },
              // Rod yield stress at room temperature
              {
                name: "f_y_rod",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("rodMaterial").mat == "cs";
                },
                desc: "Yield strength of the rod at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD_rod") - get("DI_rod")) / 2;
                  var materialType = get("rodMaterial").mat;
                  var material = get("rodMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Rod yield stress at design temperature
              {
                name: "f_yt_rod",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var rodMaterial = get("rodMaterial");
                  	return rodMaterial.mat == "cs" && rodMaterial.def == "custom";
                },
                desc: "Yield strength of the rod at design temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD_rod") - get("DI_rod")) / 2;
                  var materialType = get("rodMaterial").mat;
                  var material = get("rodMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[1];
                }
              },
              // Rod tensile strength output block
              {
                name: "TS_rod",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the rod at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("OD_rod") - get("DI_rod")) / 2;
                  var materialType = get("rodMaterial").mat;
                  var material = get("rodMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              }
            ]
            // End of dimensions menu items block
          }, 
          // End of dimensions menu block

          // Cylinder end eye menu block
          { 
            title: "Cylinder end eye",
            name: "tube_end",
            desc: "Information about the cylinder tube end fixings/end types",
            show: function show(get) {
              return get("tubeEnd") == "endEye";
            },
            items: 
            [
              // End eye figure block
              {
                name: "figureEndEye",
                type: "figure",
                src: ring,
                desc: "End eye dimensions"
              },
              // Width of tube end eye input block
              {
                name: "T_tube",
                type: "input",
                ui: "number",
                desc: "Width of tube end eye, mm",
                show: function show(get) {
                  return get("tubeEnd") == "endEye";
                }
              },
              // Outer radius of tube end eye input block
              {
                name: "R_tube",
                type: "input",
                ui: "number",
                desc: "Outer radius of tube end eye, mm",
                show: function show(get) {
                  return get("tubeEnd") == "endEye";
                }
              }, 
              // Inner diameter of tube end eye input block
              {
                name: "d_tube",
                type: "input",
                ui: "number",
                desc: "Inner diameter of tube end eye, mm",
                show: function show(get) {
                  return get("tubeEnd") == "endEye";
                }
              },
              // End eye material list input block
              {
                name: "tubeEndEyeMaterial",
                type: "input",
                desc: "Select end eye material",
                ui: uiListMaterial,
                show: function show(get) {
                  return get("tubeEnd") == "endEye";
                }
              },
              // Tube end eye 1% proof stress at room temperature
              {
                name: "R_p_tubeEndEye",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("tubeEndEyeMaterial").mat == "ss";
                },
                desc: "1% proof stress of tube end eye at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_tube");
                  var materialType = get("tubeEndEyeMaterial").mat;
                  var material = get("tubeEndEyeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Tube end eye 1% proof stress at design temperature
	      {
                name: "R_p1t_tubeEndEye",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
			var tubeEndEyeMaterial = get("tubeEndEyeMaterial");
                  	return tubeEndEyeMaterial.mat == "ss" && tubeEndEyeMaterial.def == "custom";
                },
                desc: "1% proof stress of tube end eye at design temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_tube");
                  var materialType = get("tubeEndEyeMaterial").mat;
                  var material = get("tubeEndEyeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[1];
                }
              },
	      // Tube end eye yield stress at room temperature
              {
                name: "f_y_tubeEndEye",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("tubeEndEyeMaterial").mat == "cs";
                },
                desc: "Yield strength of the tube end eye at room temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_tube");
			var materialType = get("tubeEndEyeMaterial").mat;
			var material = get("tubeEndEyeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[0];
                }
              },
	      // Tube end eye yield stress at design temperature
              {
                name: "f_yt_tubeEndEye",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var tubeEndEyeMaterial = get("tubeEndEyeMaterial");
                  	return tubeEndEyeMaterial.mat == "cs" && tubeEndEyeMaterial.def == "custom";
                },
                desc: "Yield strength of the tube end eye at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_tube");
			var materialType = get("tubeEndEyeMaterial").mat;
			var material = get("tubeEndEyeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Tensile strength of the tube end eye output block
              {
                name: "TS_tubeEndEye",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the tube end eye at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_tube");
                  var materialType = get("tubeEndEyeMaterial").mat;
                  var material = get("tubeEndEyeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // Tube end eye attachment dropdown block
              {
                name: "tubeEndEyeAttachment",
                type: "input",
                desc: "Tube end eye attachment",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "machined",
                    value: "machined"
                  }, {
                    text: "welded",
                    value: "welded"
                  }]
                },
                show: function show(get) {
                  return get("tubeEnd") == "endEye";
                }
              },
              // Tube end eye weld type dropdown block
              {
                name: "tubeEndEyeWeldType",
                type: "input",
                desc: "Tube end eye weld type",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "full-penetration",
                    value: "fPen"
                  }, {
                    text: "partial penetration",
                    value: "pPen"
                  }, {
                    text: "fillet weld",
                    value: "fWeld"
                  }]
                },
                show: function show(get) {
                  return get("tubeEndEyeAttachment") == "welded";
                }
              },
              // Fatigue choice dropdown block
              {
                name: "fatigueChoice1",
                type: "input",
                desc: "Fatigue choice",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Limit cylinder to 15 000 cycles over life span",
                    value: "fChoice1"
                  }, {
                    text: "Perform fatigue analysis",
                    value: "fChoice2"
                  }]
                },
                show: function show(get) {
                  return get("tubeEndEyeWeldType") == "pPen" || get("tubeEndEyeWeldType") == "fWeld";
                }
              },
              // Area of weld input block
              {
                name: "weld_area_tube",
                type: "input",
                ui: "number",
                desc: "Area of weld, mm^2",
                show: function show(get) {
                  return get("fatigueChoice1") == "fChoice2";
                }
              },
              // Number of cycles input block
              {
                name: "n_cycles_manufacturer",
                type: "input",
                ui: "number",
                desc: "Number of cycles",
                show: function show(get) {
                  return get("fatigueChoice1") == "fChoice2";
                }
              }
            ]
            // End of cylinder end eye menu items block
          },
          // End of cylinder end eye menu block

          // Beginning of trunnion mounted cylinder end connection menu block
          {
            title: "Trunnion mounted cylinder end connection",
            name: "trunnion",
            desc: "Information about the cylinder trunnion",
            show: function show(get) {
              return get("tubeEnd") == "trunnion";
            },
            items: 
            [
              // Tube end eye weld type dropdown block
              {
                name: "trunnionWeldType",
                type: "input",
                desc: "Tube end eye weld type",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "full-penetration",
                    value: "fPen"
                  }, {
                    text: "partial penetration",
                    value: "pPen"
                  }, {
                    text: "fillet weld",
                    value: "fWeld"
                  }]
                }
              }, 
              // Partial penetration weld figure block
              {
                name: "trunnion_weld_partial",
                type: "figure",
                src: trunnion_weld_partial,
                desc: "",
                show: function show(get) {
                  return get("trunnionWeldType") == "pPen";
                }
              },
              // Fillet weld figure block
              {
                name: "trunnion_weld_fillet",
                type: "figure",
                src: trunnion_weld_fillet,
                desc: "",
                show: function show(get) {
                  return get("trunnionWeldType") == "fWeld";
                }
              },
              // Fatigue choice for trunnion weld dropdown block
              {
                name: "fatigueChoice3",
                type: "input",
                desc: "fatigue choice for trunnion weld",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Limit cylinder to 15 000 cycles over life span",
                    value: "fChoice1"
                  }, {
                    text: "Perform fatigue analysis",
                    value: "fChoice2"
                  }]
                },
                show: function show(get) {
                  return get("trunnionWeldType") == "pPen" || get("trunnionWeldType") == "fWeld";
                }
              },
              // Throat thickness input block
              // be performed.
              {
                name: "weld_throat_trunnion",
                type: "input",
                ui: "number",
                desc: "Throat thickness, mm",
                var: "a",
                show: function show(get) {
                  return get("trunnionWeldType") == "fWeld" && get("fatigueChoice3") == "fChoice2";
                }
              }, 
              // Minimum weld leg dimension input block
              // be performed.
              {
                name: "weld_leg_trunnion",
                type: "input",
                ui: "number",
                var: "Min(L_{1},L_{2})",
                desc: "Minimum weld leg dimension, mm",
                show: function show(get) {
                  return get("trunnionWeldType") == "pPen" && get("fatigueChoice3") == "fChoice2";
                }
              },
              // Number of cycles input block
              // analysis is to be performed. 
              {
                name: "n_cycles_manufacturer_trunnion",
                type: "input",
                ui: "number",
                desc: "Number of cycles",
                show: function show(get) {
                  return get("fatigueChoice3") == "fChoice2";
                }
              }
            ]
            // End of trunnion mounted cylinder end menu items block
          },
          // End of trunnion mounted cylinder end menu block

          // Flanged cylinder end connection menu block
          {
            title: "Flanged cylinder end connection",
            name: "flange",
            desc: "Information about the cylinder flange",
            show: function show(get) {
              return (get("tubeEnd") == "flanged" || get("tubeEnd") == "flanged_located_at_end_cover");
            },
            items: 
            [
              // Flange figure block
              {
                name: "flanged_end",
                type: "figure",
                src: flanged_end,
                desc: ""
              },
              // Flange thickness input block
              itNumPos("T_flange", "The thickness of the flange, mm", "t"),
              // Number of flange bolt input block
              itNumPos("n_bolt_flange", "The number of bolts in the flange"),
              // Flange bolt circle diameter input block
              itNumPos("D_bolt_flange", "The circle diameter of the bolts in the flange, mm", "D"),
              // Flange material list block
              {
                name: "flangeMaterial",
                type: "input",
                desc: "Select flange material",
                ui: uiListMaterial
              },
              // Flange 1% proof stress at room temperature
              {
                name: "R_p_flange",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("flangeMaterial").mat == "ss";
                },
                desc: "1% proof stress of the flange at room temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_flange");
			var materialType = get("flangeMaterial").mat;
			var material = get("flangeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[0];
                }
              },
	      // Flange 1% proof stress at design temperature
              {
                name: "R_p1t_flange",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var flangeMaterial = get("flangeMaterial");
			return flangeMaterial.mat == "ss" && flangeMaterial.def == "custom";
                },
                desc: "1% proof stress of the flange at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_flange");
			var materialType = get("flangeMaterial").mat;
			var material = get("flangeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Flange yield strength at room temperature 
              {
                name: "f_y_flange",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("flangeMaterial").mat == "cs";
                },
                desc: "Yield strength of the flange at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
			var thickness = get("T_flange");
			var materialType = get("flangeMaterial").mat;
			var material = get("flangeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[0];
                }
              },
	      // Flange yield strength at design temperature 
              {
                name: "f_yt_flange",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var flangeMaterial = get("flangeMaterial");
			return flangeMaterial.mat == "cs" && flangeMaterial.def == "custom";
                },
                desc: "Yield strength of the flange at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_flange");
			var materialType = get("flangeMaterial").mat;
			var material = get("flangeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Flange tensile strength output block
              {
                name: "TS_flange",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the flange at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_flange");
                  var materialType = get("flangeMaterial").mat;
                  var material = get("flangeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // Flange weld type dropdown menu
              {
                name: "flangedWeldType",
                type: "input",
                desc: "Flanged tube end connection weld type",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "full-penetration",
                    value: "fPen"
                  }, {
                    text: "partial penetration from two sides",
                    value: "pPen"
                  }, {
                    text: "partial penetration from one side",
                    value: "pPen1"
                  }, {
                    text: "fillet weld",
                    value: "fWeld"
                  }]
                }
              },
              // Flange two sided partial weld figure block
              {
                name: "flanged_weld_partial",
                type: "figure",
                src: flange_partial,
                desc: "",
                show: function show(get) {
                  return get("flangedWeldType") == "pPen";
                }
              },
              // Flange one sided partial weld figure block
              {
                name: "flanged_weld_partial_mono",
                type: "figure",
                src: flange_partial_mono,
                desc: "",
                show: function show(get) {
                  return get("flangedWeldType") == "pPen1";
                }
              },
              // Flange fillet weld figure block
              {
                name: "flanged_weld_fillet",
                type: "figure",
                src: flange_fillet,
                desc: "",
                show: function show(get) {
                  return get("flangedWeldType") == "fWeld";
                }
              },
              // Flange fatigue choice dropdown block
              {
                name: "fatigueChoice4",
                type: "input",
                desc: "fatigue choice for flanged weld",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Limit cylinder to 15 000 cycles over life span",
                    value: "fChoice1"
                  }, {
                    text: "Perform fatigue analysis",
                    value: "fChoice2"
                  }]
                },
                show: function show(get) {
                  return get("flangedWeldType") == "pPen" || get("flangedWeldType") == "fWeld" || get("flangedWeldType") == "pPen1";
                }
              },
              // Flange throat thickness input block
              {
                name: "weld_throat_flanged",
                type: "input",
                ui: "number",
                desc: "Throat thickness, mm",
                show: function show(get) {
                  return get("flangedWeldType") == "fWeld";
                }
              },
              // Flange minimum weld leg dimension input block
              {
                name: "weld_leg_flanged",
                type: "input",
                ui: "number",
                var: "Min(L_{1},L_{2})",
                desc: "Minimum weld leg dimension, mm",
                show: function show(get) {
                  return get("flangedWeldType") == "pPen" || get("flangedWeldType") == "pPen1";
                }
              },
              // Flange number of fatigue cycles input block
              // to analysis.
              {
                name: "n_cycles_manufacturer_flanged",
                type: "input",
                ui: "number",
                desc: "Number of cycles",
                show: function show(get) {
                  return get("fatigueChoice4") == "fChoice2";
                }
              }
            ]
            // End of flanged cylinder end menu items block
          },
          // End of flanged cylinder end menu block

          // Rod end eye menu block
          {
            title: "Rod end eye",
            name: "rod_end",
            desc: "Information about the rod end fixings/end types",
            show: function show(get) {
              return get("rodEnd") == "endEye";
            },
            items: 
            [
              // Rod end eye figure block
              {
                name: "figureEndEye2",
                type: "figure",
                src: ring,
                desc: "End eye dimensions"
              },
              // Rod end eye width input block
              {
                name: "T_rod",
                type: "input",
                ui: "number",
                desc: "Width of rod end eye, mm"
              },
              // Rod end eye outer radius input block
              {
                name: "R_rod",
                type: "input",
                ui: "number",
                desc: "Outer radius of rod end eye, mm"
              },
              // Rod end eye inner radius input block
              {
                name: "d_rod",
                type: "input",
                ui: "number",
                desc: "Inner diameter of rod end eye, mm"
              },
              // Rod end eye material list block
              {
                name: "rodEndEyeMaterial",
                type: "input",
                desc: "Select end eye material",
                ui: uiListMaterial
              },
              // Rod end eye 1% proof stress at room temperature
              {
                name: "R_p_rodEndEye",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("rodEndEyeMaterial").mat == "ss";
                },
                desc: "1% proof stress of the rod end eye at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_rod");
                  var materialType = get("rodEndEyeMaterial").mat;
                  var material = get("rodEndEyeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Rod end eye 1% proof stress at design temperature
              {
                name: "R_p1t_rodEndEye",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var rodEndEyeMaterial = get("rodEndEyeMaterial");
			return rodEndEyeMaterial.mat == "ss" && rodEndEyeMaterial.def == "custom";
                },
                desc: "1% proof stress of the rod end eye at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_rod");
			var materialType = get("rodEndEyeMaterial").mat;
			var material = get("rodEndEyeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Rod end eye material yield stress at room temperature
              {
                name: "f_y_rod_end_eye",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
			return get("rodEndEyeMaterial").mat == "cs";
                },
                desc: "Yield strength of the rod end eye at room temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_rod");
			var materialType = get("rodEndEyeMaterial").mat;
			var material = get("rodEndEyeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[0];
                }
              },
	      // Rod end eye material yield stress at design temperature
              {
                name: "f_yt_rod_end_eye",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var rodEndEyeMaterial = get("rodEndEyeMaterial");
			return rodEndEyeMaterial.mat == "cs" && rodEndEyeMaterial.def == "custom";
                },
                desc: "Yield strength of the rod end eye at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("T_rod");
			var materialType = get("rodEndEyeMaterial").mat;
			var material = get("rodEndEyeMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Rod end eye material tensile strength output block
              {
                name: "TS_rod_end_eye",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the rod end eye at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("T_rod");
                  var materialType = get("rodEndEyeMaterial").mat;
                  var material = get("rodEndEyeMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // Rod end eye attachment input block
              {
                name: "rodEndEyeAttachment",
                type: "input",
                desc: "Rod end eye attachement",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "threaded",
                    value: "threaded"
                  }, {
                    text: "welded",
                    value: "welded"
                  }]
                },
                show: function show(get) {
                  return get("rodEnd") == "endEye";
                }
              },
              // Rod end eye thread diameter input block
              {
                name: "Md_rod_eye",
                type: "input",
                ui: "number",
                desc: "Rod end eye thread metric diameter, mm",
                show: function show(get) {
                  return get("rodEndEyeAttachment") == "threaded";
                }
              },
              // Rod end eye thread pitch input block
              {
                name: "xP_rod_eye",
                type: "input",
                ui: "number",
                desc: "Rod end eye thread metric pitch, mm",
                show: function show(get) {
                  return get("rodEndEyeAttachment") == "threaded";
                }
              },
              // Rod end eye thread engagement length input block
              {
                name: "Le_rod_eye",
                type: "input",
                ui: "number",
                desc: "Rod end eye thread engagement length",
                show: function show(get) {
                  return get("rodEndEyeAttachment") == "threaded";
                }
              },
              // Rod end eye weld type dropdown block
              {
                name: "rodEndEyeWeldType",
                type: "input",
                desc: "Rod end eye weld type",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "full-penetration",
                    value: "fPen"
                  }, {
                    text: "partial penetration",
                    value: "pPen"
                  }, {
                    text: "fillet weld",
                    value: "fWeld"
                  }]
                },
                show: function show(get) {
                  return get("rodEndEyeAttachment") == "welded";
                }
              },
              // Rod end eye fatigue choice dropdown block
              {
                name: "fatigueChoice2",
                type: "input",
                desc: "fatigue choice for rod",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Limit cylinder to 15 000 cycles over life span",
                    value: "fChoice1"
                  }, {
                    text: "Perform fatigue analysis",
                    value: "fChoice2"
                  }]
                },
                show: function show(get) {
                  return get("rodEndEyeWeldType") == "pPen" || get("rodEndEyeWeldType") == "fWeld";
                }
              },
              // Rod end eye area of weld input block
              {
                name: "weld_area_rod",
                type: "input",
                ui: "number",
                desc: "Area of weld, mm^2",
                show: function show(get) {
                  return get("fatigueChoice2") == "fChoice2";
                }
              },
              // Rod end eye number of fatigue cycles input block
              {
                name: "n_cycles_manufacturer_rod_eye",
                type: "input",
                ui: "number",
                desc: "Number of cycles",
                show: function show(get) {
                  return get("fatigueChoice2") == "fChoice2";
                }
              }
            ]
            // End of rod end eye menu items block
          },
          // End of rod end eye menu block

          // Threaded rod end menu block
          {
            title: "Rod end thread",
            name: "rod_end_thread",
            desc: "Information about the rod end fixings/end types",
            show: function show(get) {
              return (get("rodEnd") == "threaded_fixed" || get("rodEnd") == "threaded_pinned");
            },
            items: 
            [
              // Threaded rod end diameter input block
              {
                name: "Md_rod",
                type: "input",
                ui: "number",
                desc: "Rod thread metric diameter, mm"
              },
              // Threaded rod end pitch input block
              {
                name: "xP_rod",
                type: "input",
                ui: "number",
                desc: "Rod thread metric pitch, mm"
              },
              // Threaded rod end engagement length input block
              {
                name: "Le_rod",
                type: "input",
                ui: "number",
                desc: "Thread engagement length"
              }
            ] // End of threaded rod end menu items block
          }, // End of threaded rod end menu block

          // Piston thread menu block
          {
            title: "Piston thread",
            name: "pistonThread",
            desc: "Piston thread details",
            items: 
            [
              // Piston material list block
              {
                name: "pistonMaterial",
                type: "input",
                desc: "Select piston material",
                ui: uiListMaterial
              },
              // Piston thread diameter input block
              {
                name: "Md_piston",
                type: "input",
                ui: "number",
                desc: "Piston thread metric diameter, mm"
              },
              // Piston thread pitch input block
              {
                name: "xP_piston",
                type: "input",
                ui: "number",
                desc: "Piston thread metric pitch, mm"
              },
              // Thread engagement length input block
              {
                name: "Le_piston",
                type: "input",
                ui: "number",
                desc: "Thread engagement length, mm"
              },
              // Piston 1% proof stress at room temperature
              {
                name: "R_p_Piston",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("pistonMaterial").mat == "ss";
                },
                desc: "1% proof stress of the piston at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("Md_piston")) / 2;
                  var materialType = get("pistonMaterial").mat;
                  var material = get("pistonMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Piston 1% proof stress at design temperature
              {
                name: "R_p1t_Piston",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var pistonMaterial = get("pistonMaterial");
			return pistonMaterial.mat == "ss" && pistonMaterial.def == "custom";
                },
                desc: "1% proof stress of the piston at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = (get("DI") - get("Md_piston")) / 2;
			var materialType = get("pistonMaterial").mat;
			var material = get("pistonMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Piston yield strength at room temperature
              {
                name: "f_y_piston",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("pistonMaterial").mat == "cs";
                },
                desc: "Yield strength of the piston at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("Md_piston")) / 2;
                  var materialType = get("pistonMaterial").mat;
                  var material = get("pistonMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Piston yield strength at design temperature
              {
                name: "f_yt_piston",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var pistonMaterial = get("pistonMaterial");
			return pistonMaterial.mat == "cs" && pistonMaterial.def == "custom";
                },
                desc: "Yield strength of the piston at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = (get("DI") - get("Md_piston")) / 2;
			var materialType = get("pistonMaterial").mat;
			var material = get("pistonMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Piston tensile strength output block
              {
                name: "TS_piston",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the piston at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("Md_piston")) / 2;
                  var materialType = get("pistonMaterial").mat;
                  var material = get("pistonMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              }
            ] // End of piston menu items block
          }, // End of piston thread menu block

          // End cover menu block 
          {
            title: "End cover",
            name: "ec",
            desc: "Hydraulic cylinder end cover",
            items:
            [
              // End cover material list block
              {
                name: "ecMaterial",
                type: "input",
                desc: "Select end cover material",
                ui: uiListMaterial
              },
              // End cover attachment dropdown block
              {
                name: "ecAttachment",
                type: "input",
                desc: "End cover to cylinder tube attachment",
                ui: {
                  type: "dropdown",
                  options: [
                  /*{
                      text: "threaded",
                      value: "threaded"
                    },*/{
                    text: "Bolted",
                    value: "bolted"
                  }, {
                    text: "Welded",
                    value: "welded"
                  }]
                }
              },
              // End cover thickness input block
              itNumPos("ect", "Thickness of end cover"),
              // Steering edge length input block
              itNumPos(
                "L_steering_edge", "Length of the steering edge, mm", "L", "Look at the figure (Length of the steering edge should not exceed 3 mm)",
                function show(get) {return get("ecAttachment") == "welded";}),
              // Steering edge figure block
              {
                name: "steeringedge",
                type: "figure",
                src: steering_edge,
                desc: "Cylinder length",
                show: function show(get) {
                  return get("ecAttachment") == "welded";
                }
              }, 
              // End cover 1% proof stress at room temperature
              {
                name: "R_p_ec",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("ecMaterial").mat == "ss";
                },
                desc: "1% proof stress of the end cover at room temperature (MPa)",

                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("ect");
                  var materialType = get("ecMaterial").mat;
                  var material = get("ecMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // End cover 1% proof stress at design temperature
              {
                name: "R_pt_ec",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var ecMaterial = get("ecMaterial");
			return ecMaterial.mat == "ss" && ecMaterial.def == "custom";
                },
                desc: "1% proof stress of the end cover at design temperature (MPa)",

                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("ect");
			var materialType = get("ecMaterial").mat;
			var material = get("ecMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // End cover yield strength at room temperature
              {
                name: "f_y_ec",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("ecMaterial").mat == "cs";
                },
                desc: "Yield strength of the end cover at room temperature (MPa)",

                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("ect");
                  var materialType = get("ecMaterial").mat;
                  var material = get("ecMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // End cover yield strength at design temperature
              {
                name: "f_yt_ec",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var ecMaterial = get("ecMaterial");
			return ecMaterial.mat == "cs" && ecMaterial.def == "custom";
                },
                desc: "Yield strength of the end cover at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = get("ect");
			var materialType = get("ecMaterial").mat;
			var material = get("ecMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // End cover tensile strength output block
              {
                name: "TS_ec",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the end cover at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = get("ect");
                  var materialType = get("ecMaterial").mat;
                  var material = get("ecMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // End cover bolt material input block
              {
                name: "ecBoltMaterial",
                type: "input",
                desc: "Select end cover bolt material",
                ui: uiBoltMaterial,
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolt nominal diameter input block
              {
                name: "d_n_bolt_ec",
                type: "input",
                ui: "number",
                desc: "Bolt nominal diameter, mm",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolt thread pitch input block
              {
                name: "L_p_bolt_ec",
                type: "input",
                ui: "number",
                desc: "Bolt thread pitch, mm",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolt stress area input block
              {
                name: "A_s_bolt_ec",
                type: "input",
                ui: "number",
                desc: "Bolt stress area (per bolt), mm^2",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover number of bolts input block
              {
                name: "n_bolts_ec",
                type: "input",
                ui: "number",
                desc: "Number of bolts",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolt thread length input block
              {
                name: "l_thread_bolt_ec",
                type: "input",
                ui: "number",
                desc: "Bolt thread length, mm",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolt pitch circle diameter input block
              {
                name: "C_bolt_circle",
                type: "input",
                ui: "number",
                desc: "Bolt pitch circle diameter, mm",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                }
              },
              // End cover bolts yield strength output block
              {
                name: "R_p_ec_bolt",
                var: "R_{p1, bolt}",
                type: "output",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                },
                desc: "Yield strength of end cover bolts at room temperature (MPa)",

                ui: {
                  type: "text"
                },
                value: function value(get) {

                  return get("ecBoltMaterial").f_y;
                }

              },
              // End cover bolts tensile strength output block
              {
                name: "TS_ec_bolt",
                var: "R_{m, bolt}",
                type: "output",
                show: function show(get) {
                  return get("ecAttachment") == "bolted";
                },
                desc: "Tensile strength of the end cover bolts at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  return get("ecBoltMaterial").TS;
                }

              }
            ]
            // End of end cover menu items block
          },
          // End of end cover menu block

          // Stuffing box menu block
          {
            title: "Stuffing box",
            name: "sb",
            desc: "Stuffing box details",
            items:
            // Stuffing box menu items 
            [
              // Stuffing box material list block
              {
                name: "sbMaterial",
                type: "input",
                desc: "Select stuffing box material",
                ui: uiListMaterial
              },
              // Stuffing box 1% proof stress at room temperature
              {
                name: "R_p_sb",
                var: "R_{p1}",
                type: "output",
                show: function show(get) {
                  return get("sbMaterial").mat == "ss";
                },
                desc: "1% proof stress of the stuffing box at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("OD_rod")) / 2;
                  var materialType = get("sbMaterial").mat;
                  var material = get("sbMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Stuffing box 1% proof stress at design temperature
              {
                name: "R_p1t_sb",
                var: "R_{p1t}",
                type: "output",
                show: function show(get) {
			var sbMaterial = get("sbMaterial");
			return sbMaterial.mat == "ss" && sbMaterial.def == "custom";
                },
                desc: "1% proof stress of the stuffing box at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = (get("DI") - get("OD_rod")) / 2;
			var materialType = get("sbMaterial").mat;
			var material = get("sbMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Stuffing box yield strength at room temperature
              {
                name: "f_y_sb",
                var: "R_{eh}",
                type: "output",
                show: function show(get) {
                  return get("sbMaterial").mat == "cs";
                },
                desc: "Yield strength of the stuffing box at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("OD_rod")) / 2;
                  var materialType = get("sbMaterial").mat;
                  var material = get("sbMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[0];
                }
              },
	      // Stuffing box yield strength at design temperature
              {
                name: "f_yt_sb",
                var: "R_{et}",
                type: "output",
                show: function show(get) {
			var sbMaterial = get("sbMaterial");
			return sbMaterial.mat == "cs" && sbMaterial.def == "custom";
                },
                desc: "Yield strength of the stuffing box at design temperature (MPa)",
                ui: {
			type: "text"
                },
                value: function value(get) {
			var thickness = (get("DI") - get("OD_rod")) / 2;
			var materialType = get("sbMaterial").mat;
			var material = get("sbMaterial");
			var data = getYieldAndStrength(thickness, material);
			return data[1];
                }
              },
              // Stuffing box tensile strength output block
              {
                name: "TS_sb",
                var: "R_m",
                type: "output",
                desc: "Tensile strength of the stuffing box at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  var thickness = (get("DI") - get("OD_rod")) / 2;
                  var materialType = get("sbMaterial").mat;
                  var material = get("sbMaterial");
                  var data = getYieldAndStrength(thickness, material);
                  return data[2];
                }
              },
              // Stuffing box attachment input block
              {
                name: "sbAttachment",
                type: "input",
                desc: "Stuffing box to cylinder tube attachment",
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Threaded",
                    value: "threaded"
                  }, {
                    text: "Bolted",
                    value: "bolted"
                  }]
                }
              },
              // Stuffing box thread type input block
              {
                name: "sbThreadType",
                type: "input",
                desc: "Stuffing box thread type",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                },
                ui: {
                  type: "dropdown",
                  options: [{
                    text: "Internal",
                    value: "internal"
                  }, {
                    text: "External",
                    value: "external"
                  }]
                }
              },
              // Stuffing box figure block
              {
                name: "dufig",
                type: "figure",
                src: du,
                desc: "Cylinder length",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                }
              },
              // Stuffing box thread diameter input block
              {
                name: "Md_sb",
                type: "input",
                ui: "number",
                desc: "Stuffing box thread metric diameter, mm",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                }
              },
              // Stuffing box thread pitch input block
              {
                name: "xP_sb",
                type: "input",
                ui: "number",
                desc: "Stuffing box thread metric pitch, mm",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                }
              },
              // Stuffing box thread engagement length input block
              {
                name: "Le_sb",
                type: "input",
                ui: "number",
                desc: "Thread engagement length, mm",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                }
              },
              // Stuffing box thread relief diameter input block
              {
                name: "Du",
                type: "input",
                var: "D_u",
                ui: "number",
                desc: "Thread relief diameter, mm",
                show: function show(g) {
                  return g("sbAttachment") == "threaded";
                }
              },
              // Stuffing box bolt material list block
              {
                name: "sbBoltMaterial",
                type: "input",
                desc: "Select bolt material",
                ui: uiBoltMaterial,
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box bolt nominal diameter input block
              {
                name: "d_n_bolt",
                type: "input",
                ui: "number",
                desc: "Bolt nominal diameter, mm",
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box bolt thread input block
              {
                name: "L_p_bolt",
                type: "input",
                ui: "number",
                desc: "Bolt thread pitch, mm",
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box bolt stress area input block
              {
                name: "A_s_bolt",
                type: "input",
                ui: "number",
                desc: "Stress area bolt, mm^2",
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box number of bolts input block
              {
                name: "n_bolts",
                type: "input",
                ui: "number",
                desc: "Number of bolts",
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box threaded bolt engagement length input block
              {
                name: "l_thread_bolt",
                type: "input",
                ui: "number",
                desc: "Threaded engagement length, mm",
                show: function show(g) {
                  return g("sbAttachment") == "bolted";
                }
              },
              // Stuffing box bolt material 1% proof stress output block
              {
                name: "R_p_sb_bolt",
                var: "R_{p1, bolt}",
                type: "output",
                show: function show(get) {
                  return get("sbAttachment") == "bolted";
                },
                desc: "1% proof stress at room temperature of stuffing box bolts (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {

                  return get("sbBoltMaterial").f_y;
                }
              },
              // Stuffing box bolt material tensile strength output block
              {
                name: "TS_sb_bolt",
                var: "R_{m, bolt}",
                type: "output",
                show: function show(get) {
                  return get("sbAttachment") == "bolted";
                },
                desc: "Tensile strength of the stuffing box bolts at room temperature (MPa)",
                ui: {
                  type: "text"
                },
                value: function value(get) {
                  return get("sbBoltMaterial").TS;
                }
              }
            ]
            // End of stuffing box menu items
          }
          // End of stuffing box menu
        ]
        // End of interface menu items block
      },
      // End of interface menu block

      // Submission block
      submission: 
      {
        desc: "Please upload the design documentation as required by the DNV GL rules. \
          It will be submitted along with the output of the self-check tool to an \
          approval engineer for review.",
        docs:
        [
          //
          {
            title: "Hydraulic cylinder drawings",
            required: true,
            desc: ""
          },
          // 
          {
            title: "Design documentation",
            required: true,
            desc: "Please provide a design report (in pdf format) where you show \
              sufficient capacity of the hydraulic cylinder."
          },
          //
          {
            title: "Other relevant documents"
          }
        ]
      },
      // End of submission block

      // Rules block
      rules: 
      [
        // Non-check items rule
        {
          desc: "Warn non-checked items",
          ref: "None",
          rule: function rule(ecAttachment, sbAttachment) 
          {
            if (ecAttachment == "other") {
              this.warn([], "Please note that the end cover attachment has \
                not been checked.");
            }
            if (sbAttachment == "other") {
              this.warn([], "Please note that the stuffing box attachment \
                has not been checked.");
            }
          }
        },
        // End of non-check items rule

        // Configuration and Euler length rule
        {
          desc: "Establishes buckling configuration and euler equivalent length factors",
          ref: "CG-0194 [3.4]",
          rule: function rule(rodEnd, tubeEnd) 
          {
            var config = "";
            var eL = 0;

            if ((tubeEnd == "flanged" || tubeEnd == "flanged_located_at_end_cover") 
              && rodEnd == "threaded_fixed")
            {
              config = "fixed-fixed";
              eL = 2;
            }
            else if ((tubeEnd == "endEye" || tubeEnd == "trunnion") 
              && (rodEnd == "endEye" || rodEnd == "threaded_pinned")) 
            {
              config = "simply-simply";
              eL = 1;
            }
            else
            {
              config = "mixed";
              eL = Math.sqrt(2);
            }

            this.set("config", config);
            this.set("eL", eL);
          }
        },
        // End of configuration and Euler length rule
        
        // Assigment and calculation of variables rule
        {
          desc: "Assignment and calculation of variables",
          ref: "CG-0194 [3.4]",
          rule: function rule(DPpull, DPpush, tubeMaterial, OD, DI, rodMaterial, OD_rod, DI_rod, L1, L2max, Stroke, eL, m_cyl) {
		// General variables
		var corrosionAllowance = 0.3;
		var L = L1 + L2max;
		var L2min = L2max - Stroke;
		var L4 = L2max - Stroke;
		var Le = L / eL;
		var L1e = L1 / eL;
		var L2e_min = L2min / eL;
		var L2e_max = L2max / eL;
		var I1 = calcInertiaMoment(OD, DI);
		var I2 = calcInertiaMoment(OD_rod, DI_rod);
		var Z = calcZFactor(I1, I2, eL, L1e, L2e_max);

		var cylinder = configureCylinder(this);
    console.log("*************************************************************************")
    console.log(cylinder)

		var Pa = 0;
		if (this.get("tubeEnd") == "trunnion" && this.get("rodInside") == "yes") {
			Pa = Math.max(DPpush, DPpull) * Math.PI
				* (Math.pow(DI,2) - Math.pow(OD_rod,2)) * (0.1 / 4000);
		} else {
			Pa = DPpush * Math.PI * Math.pow(DI,2) * (0.1 / 4000);
		}
		
		// Assigning values to object variables
		this.set("cylinder", cylinder);
		this.set("Pa", Pa);
		this.set("c", corrosionAllowance);
		this.set("L", L);
		this.set("L2min", L2min);
		this.set("L4", L4);
		this.set("Le", Le);
		this.set("L1e", L1 / eL);
		this.set("L2e_min", L2min / eL);
		this.set("L2e_max", L2max / eL);
		this.set("I1", I1);
		this.set("I2", I2);
		this.set("Z", Z)
          }
        },
        // End of assignment and calculation of variables rule

        // Impact test rule
        {
          desc: "Check impact test",
          ref: "DNVGL-ST-0194 [2.1.1]",
          rule: function rule(pv_func) 
          {
            if (pv_func != "waterjet")
            {
              this.warn([], "Please remember that the materials shall be Charpy-V impact tested \
		      at minimum design temperature. Minimum specified average energy value \
		      shall be as per relevant parts of DNVGL-RU-SHIP Pt.2 Ch.2. Please refer \
		      to $ref.");
                
            }
          }
        },
        // End of impact test rule

        // Waterjet approval of manufacturer rule
        {
          desc: "Approval of manufacturer for waterjet",
          ref: "", // Add reference to DNVGL rules
          rule: function rule(pv_func)
          {
            if (pv_func == "waterjet")
            {
              this.warn([], "For hydraulic cylinders for steering gear and \
                waterjet steering approval of manufacturer is required.");
            }
          }
        },
        // End of waterjet approval of manufacturer rule

        // Waterjet working pressure rule
        {
          desc: "Check working pressure for waterjet",
          ref: "",
          rule: function rule(pv_func, wp, DPpush, DPpull) 
          {
            if (pv_func !== "waterjet")
            {
              return;
            }

            var min_required_pressure = 1.25 * wp;
            min_required_pressure = roundToDecimal(min_required_pressure, 3);

            if (min_required_pressure > DPpush) 
            {
              this.error([], "Design push pressure must be at least 1.25 times the \
                working pressure: " + min_required_pressure + " MPa > " + DPpush + " MPa.");
            }
            else
            {
              this.warn([], "Design push pressure is larger than or equal to 1.25 times \
                the working pressure:\n" + min_required_pressure + " MPa =< " + DPpush + " MPa.");
            }

            if (min_required_pressure > DPpull)
            {
              this.error([], "Design pull pressure must be at least 1.25 times the working \
                pressure: " + min_required_pressure + " MPa > " + DPpull + "MPa.");
            }
            else
            {
              this.warn([], "Design pull pressure is larger than or equal to 1.25 times \
                the working pressure: " + min_required_pressure + " MPa =< " + DPpull + " MPa.");
            }
          }
        },
        // End of waterjet working pressure rule

        // Corrosion allowance rule
        {
          desc: "Checks that the corrosion allowance is correct",
          ref: "CG-0194 [5]",
          rule: function rule(tubeMaterial, c, medium_type) 
          {
            if (medium_type == "air" && c < 1)
            {
              this.error(["c"], "The corrosion allowance for pneumatic cylinders shall \
                be at least 1 mm: " + c + "mm < " + "1 mm");
            }
            if (tubeMaterial.mat != "ss" && c < 0.3) 
            {
              this.error(["c"], "According to $ref the corrosion allowance for carbon \
                and low allow steel is usually 0.3 mm.");
            }
          }
        },
        // End of corrosion allowance rule

        // Test pressure rule
        {
          desc: "Check test pressure",
          ref: "CG-0194 sec.4 [1.1]",
          rule: function rule(p_test_pull, p_test_push, DPpush, DPpull)
          {
            var p_test_pull_req = 1.5 * DPpull;
            var p_test_push_req = 1.5 * DPpush;

            p_test_pull_req = roundToDecimal(p_test_pull_req, 3);
            p_test_push_req = roundToDecimal(p_test_push_req, 3);

            if (p_test_pull < p_test_pull_req)
            {
              this.error([], "The test pull pressure is less than 1.5 times the design pull pressure: " 
                + p_test_pull + " MPa < " + p_test_pull_req + " MPa. Please refer to $ref.");
            }
            else
            {
              this.warn([], "The test pull pressure is more than or equal to 1.5 times the design \
                pull pressure: " + p_test_pull + " MPa >= " + p_test_pull_req + " MPa. Please refer \
                to $ref.");
            }

            if (p_test_push < p_test_push_req)
            {
              this.error([], "The test push pressure is less than 1.5 times the design push pressure: " 
                + p_test_push + " MPa < " + p_test_push_req + " MPa. Please refer to $ref.");
            }
            else
            {
              this.warn([], "The test push pressure is more than or equal to 1.5 times the design \
                push pressure: " + p_test_push + " MPa >= " + p_test_push_req + " MPa. Please refer \
                to $ref.");
            }
          }
        },
        // End of test pressure rule

        // Tube sigma rule
        {
          desc: "Find sigma_tube",
          ref: "Pt.4 Ch.7 Sec.4 [2.5]", // Add DNVGL-RU-SHIP to reference?
          rule: function rule(tubeMaterial, pv_func, material_define, OD, DI) 
          {
            if (pv_func == "waterjet") {
              var factors = [1.7, 1.7, 3.5];
            } else {
              var factors = [1.7, 1.6, 2.7];
            }

            var thickness = (OD - DI) / 2;
            var materialType = tubeMaterial.mat;

            if (materialType == "E355+SR, EN10305-1") {
              thickness = OD;
            }

            var f_y_tube = tubeMaterial.f_y;
            this.set("f_y_tube", f_y_tube[0][2]);

            var sigma_tube = calc_nominal_stress(thickness, factors, tubeMaterial);
            this.set("sigma_tube", sigma_tube);
          }
        },
        // End of tube sigma rule

        // Rod sigma rule
        {
          desc: "Find sigma_rod",
          ref: "Pt.4 Ch.7 Sec.4 [2.5]",
          rule: function rule(rodMaterial, OD_rod, DI_rod, material_define) 
          {
            var thickness = (OD_rod - DI_rod) / 2;
            var factors = [1.7, 1.6, 2.7];
            var materialType = rodMaterial.mat;
            var sigma_rod = calc_nominal_stress(thickness, factors, rodMaterial);
            this.set("sigma_rod", sigma_rod);
          }
        },
        // End of rod sigma rule

        // End cover sigma rule
        {
          desc: "Find sigma_ec",
          ref: "Pt.4 Ch.7 Sec.4 [2.5]",
          rule: function rule(ecMaterial, pv_func, ect, material_define) 
          {
            var thickness = ect;
            if (pv_func == "waterjet") {
              var factors = [1.7, 1.7, 3.5];
            } else {
              var factors = [1.7, 1.6, 2.7];
            }
            var materialType = ecMaterial.mat;
            var sigma_ec = calc_nominal_stress(thickness, factors, ecMaterial);
            this.set("sigma_ec", sigma_ec);
          }
        },
        // End of end cover sigma rule

        // Required minimum tube thickness rule
        {
          desc: "Calculates required minimum tube thickness",
          ref: "CG-0194 [5]",
          "rule": function rule(DPpush, DPpull, sigma_tube, v, DI, OD, c) 
          {
            var DP_max = Math.max(DPpush, DPpull);
            this.set("DP_max", DP_max);
            var s_nom = (OD - DI) / 2;
            this.set("s_nom", s_nom);
            var s_min = DP_max * OD / (20 * sigma_tube * v + DP_max) + c;
            this.set("s_min", s_min);

            s_min = roundToDecimal(s_min, 3);
            s_nom = roundToDecimal(s_nom, 3);
            if (!(s_nom >= s_min)) {
              this.error(["s_nom"], "According to $ref nominal tube thickness is not sufficient: " 
                + s_min + " mm > " + s_nom + " mm.");
            } else {
              this.warn(["s_nom"], "According to $ref nominal tube thickness is sufficient: "
                + s_min + " mm =< " + s_nom + " mm.");
            }
          }
        },
        // End of required minimum tube thickness rule

        // Buckling safety factor rule
        {
		desc: "Calculates buckling safety factor [simple]",
		ref: "DNVGL-ST-0194 [3.2]",
		"rule": function rule(cylinder) {
			var requiredSafetyFactor = 4.0;
			var safetyFactor = undefined;
			var stroke = cylinder.piston.stroke;
			var data = undefined;
			var prefix = "The buckling safety factor is: ";
			var condition = "";
			var reference = "Please refer to $ref.";
			var validCurve = true;
			var curveError = "Insignificant input to perform buckling calculations. \
				Force-stroke and pressure-stroke curves need at least two entries.";
			if (cylinder.bucklingMethod == "simple_static") {
				safetyFactor = bucklingStatic(cylinder);
			} else if (cylinder.bucklingMethod == "simple_force") {
				data = bucklingForceCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
				validCurve = (typeof safetyFactor != "undefined");
			} else if (cylinder.bucklingMethod == "simple_pressure") {
				data = bucklingPressureCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
				validCurve = (typeof safetyFactor != "undefined");
			} else {
				return;
			}
			if (safetyFactor < requiredSafetyFactor) {
				condition = safetyFactor + " < " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.error([], prefix + condition + reference);
			} else if (safetyFactor >= requiredSafetyFactor) {
				condition = safetyFactor + " >= " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.warn([], prefix + condition + reference);
			} else if (!validCurve) {
				this.error([], curveError);
			}
	  } 
        },
        // End of buckling safety factor rule

        // Accurate buckling safety factor rule
	{
		desc: "Calculates buckling safety factor [accurate]",
		ref: "DNVGL-ST-0194 [A.4]",
		"rule": function rule(cylinder) {
            		var requiredSafetyFactor = 2.7;
			var safetyFactor = undefined;
			var stroke = cylinder.piston.stroke;
			var data = undefined;
			var prefix = "The buckling safety factor is: ";
			var condition = "";
			var reference = "Please refer to $ref.";
			var validCurve = true;
			var curveError = "Insignificant input to perform buckling calculations. \
				Force-stroke and pressure-stroke curves need at least two entries.";
			if (cylinder.bucklingMethod == "accurate_static") {
				safetyFactor = bucklingStatic(cylinder);
			} else if (cylinder.bucklingMethod == "accurate_force") {
				data = bucklingForceCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
				validCurve = (typeof safetyFactor != "undefined");
			} else if (cylinder.bucklingMethod == "accurate_pressure") {
				data = bucklingPressureCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
				validCurve = (typeof safetyFactor != "undefined");
			} else {
				return;
			}
			if (safetyFactor < requiredSafetyFactor) {
				condition = safetyFactor + " < " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.error([], prefix + condition + reference);
			} else if (safetyFactor >= requiredSafetyFactor) {
				condition = safetyFactor + " >= " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.warn([], prefix + condition + reference);
			} else if (!validCurve) {
				this.error([], curveError);
			}
		}
        },
        // End of accurate buckling safety factor rule
	
	// EN buckling safety factor rule
	{
		desc: "Calculates buckling safety factor [EN]",
          	ref: "DNVGL-ST-0194 [A.5]",
		"rule": function rule(cylinder) {
			var requiredSafetyFactor = 2.0;
			var safetyFactor = undefined;
			var stroke = cylinder.piston.stroke;
			var data = undefined;
			var prefix = "The buckling safety factor is: ";
			var condition = "";
			var reference = "Please refer to $ref.";
			var validCurve = true;
			var curveError = "Insignificant input to perform buckling calculations. \
				Force-stroke and pressure-stroke curves need at least two entries.";
			if (cylinder.bucklingMethod == "en_static") {
				safetyFactor = bucklingStatic(cylinder);
			} else if (cylinder.bucklingMethod == "en_force") {
				data = bucklingForceCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
				validCurve = (typeof safetyFactor != "undefined");
			} else if (cylinder.bucklingMethod == "en_pressure") {
				data = bucklingPressureCurve(cylinder);
				safetyFactor = data.safetyFactor;
				stroke = data.stroke;
			} else {
				return;
			}
			if (safetyFactor < requiredSafetyFactor) {
				condition = safetyFactor + " < " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.error([], prefix + condition + reference);
			} else if (safetyFactor >= requiredSafetyFactor) {
				condition = safetyFactor + " >= " + requiredSafetyFactor + " at stroke "
					+ stroke + " mm. ";
				this.warn([], prefix + condition + reference);
			} else if (typeof safetyFactor == "undefined") {
				this.error([], curveError);
			}
		}
	}, // End of EN buckling safety factor rule

        // Rod end eye stress rule
        {
          desc: "Calculates rod end eye stress",
          ref: "CG-0194 [8]",
          "rule": function rule(T_rod, R_rod, d_rod, rodEndEyeMaterial, DPpush, DPpull, 
            DI, OD_rod, rodEnd, material_define) 
          {
            var F_pull = DPpull * 0.1 * Math.PI * (DI * DI - OD_rod * OD_rod) / 4; // unit: N
            this.set("F_pull", F_pull);

            if (rodEnd == "endEye") {
              var D_rod = R_rod * 2;
              
              //var F_push = DPpush * 0.1 * Math.PI * DI * DI / 4; // unit: N
              //this.set("F_push", F_push)
              //var F_max = Math.max(F_push, F_pull)
              //this.set("F_max", F_max)
              var t_rod = F_pull * Math.sqrt(D_rod / d_rod * (D_rod / d_rod) 
                - D_rod / d_rod + 1) / (T_rod * (D_rod - d_rod));
              var thickness = T_rod;
              var material = getYieldAndStrength(thickness, rodEndEyeMaterial);
              var sigma_rod_end_eye = material[0];
              this.set("sigma_rod_end_eye", sigma_rod_end_eye);

              sigma_rod_end_eye = roundToDecimal(sigma_rod_end_eye, 3);
              t_rod = roundToDecimal(t_rod, 3);
              this.set("t_rod", t_rod);

              if (!(t_rod < sigma_rod_end_eye)) {
                this.error(["t_rod"], "The stress in the rod end eye is above the nominal design \
                  stress: " + t_rod + " MPa > " + sigma_rod_end_eye + " MPa. Please refer to $ref.");
              } else {
                this.warn(["t_rod"], "The stress in the rod end eye is lower than or equal or \
                  equal to the nominal design stress: " + t_rod + " MPa =< " + sigma_rod_end_eye 
                  + " MPa. Please refer to $ref.");
              }
            }
          }           
        },
        // End of rod end eye stress rule

        // Tube end eye stress rule
        {
          desc: "Calculates tube end eye stress",
          ref: "CG-0194 [8]",
          "rule": function rule(T_tube, R_tube, d_tube, tubeEndEyeMaterial, DPpush, DPpull, 
            DI, OD_rod, tubeEnd, F_pull, material_define) 
          {
            var D_tube = R_tube * 2;
            var t_tube = F_pull * Math.sqrt(D_tube / d_tube * (D_tube / d_tube) 
              - D_tube / d_tube + 1) / (T_tube * (D_tube - d_tube));
            t_tube = roundToDecimal(t_tube, 3); // (1000 * t_tube) /1000;
            this.set("t_tube", t_tube);
            
            if (tubeEnd == "endEye") {

              var thickness = T_tube;
              var material = getYieldAndStrength(thickness, tubeEndEyeMaterial);
              var sigma_tube_end_eye = material[0];
              this.set("sigma_tube_end_eye", sigma_tube_end_eye);

              if (!(t_tube < sigma_tube_end_eye)) {
                this.error(["t_tube"], "The stress in the tube end eye is above the nominal \
                  design stress: " + t_tube + " MPa > " + sigma_tube_end_eye 
                  + " MPa. Please refer to $ref.");
              } else {
                this.warn(["t_tube"], "The stress in the tube end eye is lower than or \
                  equal the nominal design stress: " + t_tube + " MPa =< " + sigma_tube_end_eye 
                  + " MPa. Please refer to $ref.");
              }
            }
          }
        },
        // End of tube end eye stress rule

        // Piston thread rule
        {
          desc: "piston thread calculation",
          ref: "CG-0194 [7] piston thread",
          rule: function rule(Md_piston, xP_piston, Le_piston, DI, OD_rod, DI_rod, 
            pistonMaterial, rodMaterial, DPpush, DPpull) 
          {
            var thickness1 = (DI - Md_piston) / 2;
            var thickness2 = (OD_rod - DI_rod) / 2;
            var material1 = getYieldAndStrength(thickness1, pistonMaterial);
            var material2 = getYieldAndStrength(thickness2, rodMaterial);

            var p_design = Math.max(DPpush, DPpull);

            rule_thread_stress(this.set, this.error, this.warn, "piston", Md_piston, xP_piston, 
              Le_piston, DI, OD_rod, material1, material2, p_design, 1, true, 1);
          }
        },
        // End of piston thread rule

        // Rod end eye thread rule
        {
          desc: "Rod end eye thread calculation",
          ref: "CG-0194 [7] rod end eye thread",
          rule: function rule(rodEndEyeAttachment, Md_rod_eye, xP_rod_eye, Le_rod_eye, DI, 
            OD_rod, DI_rod, T_rod, rodEndEyeMaterial, rodMaterial, F_pull, rodEnd, DPpush, 
            DPpull, material_define)
          {  
            if (rodEnd == "endEye" && rodEndEyeAttachment == "threaded") {
              var thickness1 = T_rod;
              var thickness2 = (OD_rod - DI_rod) / 2;
              var material1 = getYieldAndStrength(thickness1, rodEndEyeMaterial);
              var material2 = getYieldAndStrength(thickness2, rodEndEyeMaterial);

              rule_thread_stress(this.set, this.error, this.warn,"rod_end_eye_push", Md_rod_eye, xP_rod_eye, 
                Le_rod_eye, DI, OD_rod, material1, material2, DPpush, 0, 0, 0);
              rule_thread_stress(this.set, this.error, this.warn, "rod_end_eye_pull", Md_rod_eye, xP_rod_eye, 
                Le_rod_eye, DI, OD_rod, material1, material2, DPpull, 0, 0, 1);
            }
          }
        },
        // End of rod end eye thread rule

        // Rod end thread rule
        {
          desc: "Rod end thread calculation",
          ref: "CG-0194 [7] rod end thread",
          rule: function rule(Md_rod, xP_rod, Le_rod, DI, OD_rod, DI_rod, pistonMaterial, 
            rodMaterial, F_pull, rodEnd, DPpush, DPpull, material_define) 
          {
            
            if (rodEnd == "threaded_fixed" || rodEnd == "threaded_pinned") {
              var thickness1 = (OD_rod - DI_rod) / 2;
              var material1 = getYieldAndStrength(thickness1, rodMaterial);

              rule_thread_stress(this.set, this.error, this.warn, "rod_thread_push", Md_rod, xP_rod, 
                Le_rod, DI, OD_rod, material1, material1, DPpush, 0, 0, 0);
              rule_thread_stress(this.set, this.error, this.warn, "rod_thread_pull", Md_rod, xP_rod, 
                Le_rod, DI, OD_rod, material1, material1, DPpull, 0, 0, 1);
            }
          }
        },
        // End of rod end thread rule

        // Stuffing box thread rule
        {
          desc: "stuffing box thread calculation",
          ref: "CG-0194 [7] stuffing box thread",
          rule: function rule(sbThreadType, Du, Le_sb, OD, DI_rod, Md_sb, xP_sb, DI, OD_rod, 
            sbMaterial, tubeMaterial, F_pull, sbAttachment, DPpush, DPpull, material_define)
          {
            var _this = this;

            if (sbAttachment == "threaded") {

              var thickness1 = (OD - DI) / 2;
              var thickness2 = (DI - OD_rod) / 2;
              var material1 = getYieldAndStrength(thickness1, sbMaterial);
              var material2 = getYieldAndStrength(thickness2, tubeMaterial);

              rule_thread_stress(this.set, this.error, this.warn, "stuffing box", Md_sb, xP_sb, Le_sb, 
                DI, OD_rod, material1, material2, DPpull, 0, 0, 1);
              // Check the flaring stress using EN14359 5.6.1.6
              // We do two checks, for push and pull pressure when for the push pressure 
              // we use the area Di^2 => Di^2 - di^2

              var De = OD,
                  Le = Le_sb,
                  Dp = Md_sb - 0.6495 * xP_sb,
                  factors = [1.5, 1.5, 2.4],
                  f = Math.min(material2[0] / 1.5, material2[2] / 2.4);
              f = roundToDecimal(f, 3);
              var ok = true;
              this.set("f_flaring", f);
              this.set("Dp", Dp);

              var func = function func(error, warn, name, type, PS, Di, Di2) {
                var v = 0.3;

                var Au = void 0,
                    Dm = void 0,
                    Et = void 0;
                if (type == "internal") {
                  Au = Math.PI * (De * De - Du * Du) / 4;
                  Et = (De - Du) / 2;
                  Dm = (De + Du) / 2;
                  var thread_type = "internal";
                } else {
                  Au = Math.PI * (Du * Du - Di * Di) / 4;
                  Et = (Du - Di) / 2;
                  Dm = (Di + Du) / 2;
                  var thread_type = "external";
                }

                _this.set("Au_" + name, Au);
                _this.set("Et_" + name, Et);
                _this.set("Dm_" + name, Dm);

                var beta = Math.pow(12 * (1 - v * v) / (Dm * Dm * Et * Et), 0.25);
                var alpha = 1 / (4 * beta * Le) * (1 + 4 * Math.exp(-beta * Le) * Math.sin(beta * Le) 
                  - Math.exp(-2 * beta * Le) * (Math.sin(2 * beta * Le) + Math.cos(2 * beta * Le)));

                _this.set("beta_" + name, beta);
                _this.set("alpha_" + name, alpha);
                var Fe = Di2 * PS * Math.abs(Dm - Dp) / (8 * Dm) / 10;
                var fb = Math.PI * Di2 * PS / 4 / Au / 10 + 6 * Fe * alpha / Math.pow(Et, 2);

                _this.set("Fe_" + name, Fe);
                _this.set("fb_" + name, fb);
                fb = roundToDecimal(fb, 3);
                if (fb > f) ok = false;
                if (!ok) {
                  error([], "The flaring stress in the " + type + " thread for the stuffing \
                    box is higher than the allowable flaring stress: " + fb + " MPa > " + f 
                    + " MPa. Please refer to EN 14359:2006 5.6.1.6.");
                } else {
                  warn([], "The flaring stress in the " + type + " thread for the stuffing \
                    box is lower than or equal to the allowable flaring stress: " + fb 
                    + " MPa =< " + f + " MPa. Please refer to EN 14359:2006 5.6.1.6.");
                }
              };
              
              func(this.error, this.warn, "pull", sbThreadType, DPpull, DI, DI * DI - OD_rod * OD_rod);
            }
          }
        },
        // End of stuffing box thread rule

        // End cover thickness rule
        {
          desc: "End cover thickness calculation",
          ref: "EN 13445-3", //"CG-0194 [6]",
          rule: function rule(ecAttachment, DPpush, sigma_ec, sigma_tube, s_nom, DI, 
            ect, OD, tubeEnd)
          {
            if (ecAttachment != "welded") 
            {
              return;
            }
            
            var fmin = Math.min(sigma_ec, sigma_tube);
            this.set("fmin", fmin);
            
            var P = DPpush * 0.1;
            var B1 = 1 - 3 * fmin * Math.pow(s_nom / (DI + s_nom), 2) / P 
              + 3 / 16 * Math.pow(DI / (DI + s_nom), 4) * P / fmin 
              - 3 / 4 * ((2 * DI + s_nom) * s_nom * s_nom) / Math.pow(DI + s_nom, 3);
            var A1 = B1 * (1 - B1 * s_nom / (2 * (DI + s_nom)));
            var C1a = 0.40825 * A1 * (DI + s_nom) / DI;
            var C1b = 0.299 * (1 + 1.7 * s_nom / DI);
            var C1 = Math.max(C1a, C1b);
            var C2 = calc_C2(DI, s_nom, 0.3, fmin, DPpush);

            this.set("A1", A1);
            this.set("B1", B1);
            this.set("C1", C1);
            this.set("C2", C2);

            var e_ec1 = C1 * DI * Math.sqrt(P / sigma_ec);
            var e_ec2 = C2 * DI * Math.sqrt(P / fmin);
            var e_ec = e_ec1;

            if (C2 >= 0.3 && e_ec1 < e_ec2) 
            {
              e_ec = e_ec2;
            }

            e_ec = roundToDecimal(e_ec, 3);
            this.set("e_ec", e_ec);
            
            if (ect < e_ec) {
              this.error([], "The end cover thickness is not sufficient: " 
                + ect + " mm < " + e_ec + " mm. Please refer to $ref.");
            } else {
              this.warn([], "The end cover thickness is sufficient: " 
                + ect + " mm >= " + e_ec + " mm. Please refer to $ref.");
            }

            var ratio_do_di = OD / DI;
            this.set("ratio_do_di", ratio_do_di);

            var sd_end_cover = (DPpush/10) * Math.pow(((341 / 350) - (3 / 7) * ((DI + 2 * s_nom) / DI)) * (DI / ect), 2)
            this.set("sd_end_cover", sd_end_cover);
          }
        },
        // End cover thickness rule

        // Tube end welding rule
        {
          desc: "Welding check tube end eye",
          ref: "CG-0194 [10] & RP-C203",
          rule: function rule(weld_area_tube, F_pull, Pa, fatigueChoice1, tubeEndEyeWeldType, 
            n_cycles_manufacturer, tubeEndEyeAttachment)
          {
            if (tubeEndEyeAttachment == "welded")
            {
              calc_fatigue(this.set, this.error, this.warn, "tube_end_eye", Pa, F_pull, tubeEndEyeWeldType, 
                fatigueChoice1, weld_area_tube, n_cycles_manufacturer);
            }
          }
        },
        // End of tube end welding rule

        // Rod end eye welding rule
        {
          desc: "Welding check rod end eye",
          ref: "CG-0194 [10] & RP-C203",
          rule: function rule(weld_area_rod, F_pull, Pa, fatigueChoice2, rodEndEyeWeldType, 
            n_cycles_manufacturer_rod_eye, rodEndEyeAttachment)
          {
            if (rodEndEyeAttachment == "welded")
            {
              calc_fatigue(this.set, this.error, this.warn, "rod_end_eye", Pa, F_pull, rodEndEyeWeldType, 
                fatigueChoice2, weld_area_rod, n_cycles_manufacturer_rod_eye);
            }
          }
        },
        // End of rod end eye welding rule

        // Trunnion welding rule
        {
          desc: "Welding check trunnion",
          ref: "CG-0194 [10] & RP-C203",
          rule: function rule(weld_throat_trunnion, weld_leg_trunnion, F_pull, Pa, fatigueChoice3, 
            trunnionWeldType, n_cycles_manufacturer_trunnion, tubeEnd, OD) 
          {
            if (tubeEnd == "trunnion")
            {
              if (trunnionWeldType == "fWeld")
              {
                var weld_area_trunnion = 2 * weld_throat_trunnion * OD * Math.PI;
              }
              else
              {
                var weld_area_trunnion = 2 * weld_leg_trunnion * OD * Math.PI;
                this.set("weld_area_trunnion", weld_area_trunnion);
              }
              calc_fatigue(this.set, this.error, this.warn, "trunnion", Pa, F_pull, trunnionWeldType, 
                fatigueChoice3, weld_area_trunnion, n_cycles_manufacturer_trunnion);
            }
          }
        },
        // End of trunnion welding rule

        // Stuffing box bolt rule
        {
          desc: "Bolted stuffing box check",
          ref: "EN 13445-3 & EN 14359",
          rule: function rule(d_n_bolt, L_p_bolt, sbBoltMaterial, tubeMaterial, A_s_bolt, n_bolts, 
            l_thread_bolt, sbAttachment, F_pull, DI, OD, tubeEnd, Pa)
          {
            if (sbAttachment == "bolted")
            {
              var thickness = (OD - DI) / 2;
              var materialType = tubeMaterial.mat;
              if (materialType == "E355+SR, EN10305-1")
              {
                thickness = OD;
              }
              var material = getYieldAndStrength(thickness, tubeMaterial);
              this.set("material", material);

              calc_bolts(this.set, this.error, this.warn, "stuffing_box", d_n_bolt, L_p_bolt, sbBoltMaterial,
                material, A_s_bolt, n_bolts, l_thread_bolt, sbAttachment, F_pull, tubeEnd, Pa);
            }
          }
        },
        // End of stuffing box bolt rule

        // End cover bolt length rule
        {
          desc: "Check end cover bolt length",
          ref: "EN 13445-3 [10.5.3]",
          rule: function rule(d_n_bolt_ec, L_p_bolt_ec, ecBoltMaterial, A_s_bolt_ec, n_bolts_ec, 
            l_thread_bolt_ec, ecAttachment, DPpush, C_bolt_circle, sigma_ec, ect, F_pull, tubeMaterial, 
            OD, DI, tubeEnd, Pa)
          {
            if (ecAttachment == "bolted")
            {
              var thickness = (OD - DI) / 2;
              var materialType = tubeMaterial.mat;

              if (materialType == "E355+SR, EN10305-1")
              {
                thickness = OD;
              }
              
              var material = getYieldAndStrength(thickness, tubeMaterial);

              calc_bolts(this.set, this.error, this.warn, "end_cover", d_n_bolt_ec, L_p_bolt_ec, ecBoltMaterial, 
                material, A_s_bolt_ec, n_bolts_ec, l_thread_bolt_ec, ecAttachment, F_pull, tubeEnd, Pa);

              var e_req_end_cover_bolted = 0.41 * C_bolt_circle * Math.sqrt(DPpush / 10 / sigma_ec);
              e_req_end_cover_bolted = roundToDecimal(e_req_end_cover_bolted, 3);
              this.set("e_req_end_cover_bolted", e_req_end_cover_bolted);

              if (e_req_end_cover_bolted > ect)
              {
                this.error([], "The thickness of the bolted end cover is lower than the required thickness: " 
                  + ect + " mm < " + e_req_end_cover_bolted + " mm. Please refer to $ref.");
              }
              else
              {
                this.warn([], "The thickness of the bolted end cover is higher than or equal to the required \
                  thickness: " + ect + " mm >= " + e_req_end_cover_bolted + " mm. Please refer to $ref.");
              }
            }
          }
        }, // End of end cover bolt length rule

        // Steering edge length rule
        {
          desc: "Check steering edge length",
          ref: "CG-0194 sec.2 [9.2]",
          rule: function rule(L_steering_edge, ecAttachment) 
          {
            if (ecAttachment == "welded")
            {
              if (L_steering_edge > 3)
              {
                this.error([], "The steering edge between the end cover and tube shall not be greater \
                  than 3 mm, when hydraulic cylinder is designed for more than 15000 full load cycles: " 
                  + L_steering_edge + " mm > 3 mm." + " Please refer to $ref.");
              }
              else 
              {
                this.warn([], "The steering edge between the end cover and tube is not greater than 3 mm, \
                  when hydraulic cylinder is designed for more than 15000 full load cycles: " 
                  + L_steering_edge + " mm =< 3 mm." + " Please refer to $ref.");
              }
            }
          }
        },
        // End of steering edge length rule

        // Full penetration flange weld rule
        {
          desc: "Check flange weld for full penetration",
          ref: "DNVGL-CG-0194",
          rule: function rule(tubeEnd, flangeMaterial, T_flange, n_bolt_flange, D_bolt_flange, Pa, OD, 
            flangedWeldType, weld_leg_flanged, F_pull)
          {
            if ((tubeEnd == "flanged" || tubeEnd == "flanged_located_at_end_cover") && flangedWeldType == "fPen")
            {
              var thickness = T_flange;
              var materialType = flangeMaterial;
              var material = flangeMaterial;
              var data = getYieldAndStrength(thickness, material);

              var sigma_flange = data[0] / 1.5;
              sigma_flange = roundToDecimal(sigma_flange, 3);
              this.set("sigma_flange", sigma_flange);
              var L_moment = (D_bolt_flange - OD) / 2;
              var F_per_segment = Math.max(Pa,F_pull*0.001) / n_bolt_flange * 1000;

              var M_flange = F_per_segment * L_moment;
              var S_cyl = Math.PI * OD / n_bolt_flange;
              var W_flange_segment = S_cyl * Math.pow(T_flange, 2) / 6;
              var A_flange_segment = S_cyl * T_flange;

              var sigma_bending = M_flange / W_flange_segment;
              this.set("sigma_bending", sigma_bending);

              var tau_shear = F_per_segment / A_flange_segment;
              this.set("tau_shear", tau_shear);

              var sigma_eqi = Math.sqrt(Math.pow(sigma_bending, 2) + 3 * Math.pow(tau_shear, 2));
              sigma_eqi = roundToDecimal(sigma_eqi, 3);
              this.set("sigma_eqi", sigma_eqi);

              if (sigma_eqi > sigma_flange)
              {
                this.error([], "The equivalent stress in the flange is higher than the allowable \
                  stress: " + sigma_eqi + " MPa > " + sigma_flange + " MPa. Please refer to $ref.");
              }
              else
              {
                this.warn([], "The equivalent stress in the flange is lower than or equal the \
                allowable stress: " + sigma_eqi + " MPa  =< " + sigma_flange + " MPa. Please \
                refer to $ref.");
              }
            }
          }
        },
        // End of full penetration flange weld rule

        // Fillet flange weld rule
        {
          desc: "Check flange weld for fillet",
          ref: "DNVGL-CG-0194",
          rule: function rule(DPpush, OD, DI, Pa, F_pull, flangeMaterial, tubeEnd, weld_throat_flanged, 
            flangedWeldType, T_flange, n_bolt_flange, D_bolt_flange, weld_leg_flanged, n_cycles_manufacturer_flanged)
          {
            if ((tubeEnd == "flanged" || tubeEnd == "flanged_located_at_end_cover") && flangedWeldType == "fWeld")
            {
              var thickness = T_flange;
              var materialType = flangeMaterial;
              var material = flangeMaterial;
              var data = getYieldAndStrength(thickness, material);
              var sigma_flange = Math.min(data[0] / 1.7, data[2] / 2.7);

              sigma_flange = roundToDecimal(sigma_flange, 3);
              this.set("sigma_flange", sigma_flange);

              var L_moment = (D_bolt_flange - OD) / 2;
              var F_per_segment = Pa / n_bolt_flange * 1000;
              var M_flange = F_per_segment * L_moment;
              var S_cyl = Math.PI * OD / n_bolt_flange;
              var I_flange_segment = S_cyl * (Math.pow(weld_throat_flanged, 3) / 12 
                + 2 * weld_throat_flanged * Math.pow(weld_throat_flanged / 2 + T_flange / 2, 2));
              var W_flange_segment = 2 * I_flange_segment / (weld_throat_flanged + T_flange);
              var A_flange_segment = 2 * S_cyl * weld_throat_flanged;

              var sigma_bending = M_flange / (Math.sqrt(2) * W_flange_segment);
              this.set("sigma_bending", sigma_bending);

              var tau_bending = sigma_bending;
              this.set("tau_bending", tau_bending);

              var tau_shear = F_per_segment / A_flange_segment;
              this.set("tau_shear", tau_shear);

              var sigma_eqi = Math.sqrt(Math.pow(sigma_bending, 2) + 3 * Math.pow(tau_shear, 2) 
                + 3 * Math.pow(tau_bending, 2));
              sigma_eqi = roundToDecimal(sigma_eqi, 3);
              this.set("sigma_eqi", sigma_eqi);

              if (sigma_eqi > sigma_flange)
              {
                this.error([], "The equivalent stress in the flange from pushing force is higher \
                  than the allowable stress: " + sigma_eqi + " MPa > " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }
              else
              {
                this.warn([], "The equivalent stress in the flange from pushing force is lower than \
                  or equal the allowable stress: \n" + sigma_eqi + " MPa =< " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }

              var F_per_segment_2 = F_pull / n_bolt_flange;
              var M_flange_2 = F_per_segment_2 * L_moment;

              var sigma_bending_2 = M_flange_2 / (Math.sqrt(2) * W_flange_segment);
              this.set("sigma_bending_2", sigma_bending_2);

              var tau_bending_2 = sigma_bending;
              this.set("tau_bending_2", tau_bending_2);

              var tau_shear_2 = F_per_segment_2 / A_flange_segment;
              this.set("tau_shear_2", tau_shear_2);

              var sigma_eqi_2 = Math.sqrt(Math.pow(sigma_bending_2, 2) 
                + 3 * Math.pow(tau_shear_2, 2) + 3 * Math.pow(tau_bending_2, 2));
              sigma_eqi_2 = roundToDecimal(sigma_eqi_2, 3);
              this.set("sigma_eqi_2", sigma_eqi_2);
              
              if (sigma_eqi_2 > sigma_flange)
              {
                this.error([], "The equivalent stress in the flange from pulling force is higher than \
                  the allowable stress: " + sigma_eqi_2 + " MPa > " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }
              else
              {
                this.warn([], "The equivalent stress in the flange from pulling force is lower than or \
                  equal the allowable stress: " + sigma_eqi_2 + " MPa =< " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }

              var sigma_eqi_fatigue_push = Math.sqrt(Math.pow(sigma_bending, 2) 
                + 0.2 * Math.pow(tau_shear, 2) + Math.pow(tau_bending, 2));
              this.set("sigma_eqi_fatigue_push", sigma_eqi_fatigue_push);

              var sigma_eqi_fatigue_pull = Math.sqrt(Math.pow(sigma_bending_2, 2) 
                + 0.2 * Math.pow(tau_shear_2, 2) + Math.pow(tau_bending_2, 2));
              this.set("sigma_eqi_fatigue_pull", sigma_eqi_fatigue_pull);

              var delta_stress = sigma_eqi_fatigue_push + sigma_eqi_fatigue_pull;
              this.set("delta_stress", delta_stress);

              var m_2 = 3;
              var log_a = 10.97;
              var log_N = log_a - m_2 * Math.log10(delta_stress);
              var n_fatigue = Math.pow(10, log_N);

              n_fatigue = Math.round(n_fatigue);
              
              if (n_fatigue < n_cycles_manufacturer_flanged)
              {
                this.error([], "The flange's fillet weld weld has a predicted number of cycles to \
                  failure of: " + n_fatigue.toFixed(0) + " < " + n_cycles_manufacturer_flanged + ". \
                  Please refer to DNVGL-RP-C203");
              }
              else
              {
                this.warn([], "The flange's fillet weld has a predicted number of cycles to \
                  failure of: " + n_fatigue.toFixed(0) + " >= " + n_cycles_manufacturer_flanged + ". \
                  Please refer to DNVGL-RP-C203");
              }
            }
          }
        }, // End of fillet flange weld rule

        // Partial penetration flange weld rule
        {
          desc: "Check flange weld for partial penetration",
          ref: "DNVGL-CG-0194",
          rule: function rule(DPpush, OD, DI, Pa, F_pull, flangeMaterial, tubeEnd, flangedWeldType, T_flange, 
            n_bolt_flange, D_bolt_flange, weld_leg_flanged, n_cycles_manufacturer_flanged)
          {
            if ((tubeEnd == "flanged" || tubeEnd == "flanged_located_at_end_cover") 
              && (flangedWeldType == "pPen" || flangedWeldType == "pPen1"))
            {
              var thickness = T_flange;
              var materialType = flangeMaterial;
              var material = flangeMaterial;
              var data = getYieldAndStrength(thickness, material);

              var sigma_flange = Math.min(data[0] / 1.7, data[2] / 2.7);
              sigma_flange = roundToDecimal(sigma_flange, 3);
              this.set("sigma_flange", sigma_flange);

              if (flangedWeldType == "pPen") 
              {
                T_flange = T_flange - 2 * weld_leg_flanged;
              }
              else
              {
                T_flange = T_flange - weld_leg_flanged;
              }

              var L_moment = (D_bolt_flange - OD) / 2;
              var F_per_segment = Pa / n_bolt_flange * 1000;
              var M_flange = F_per_segment * L_moment;
              var S_cyl = Math.PI * OD / n_bolt_flange;
              var I_flange_segment = S_cyl * (Math.pow(weld_leg_flanged, 3) / 12 
                + 2 * weld_leg_flanged * Math.pow(weld_leg_flanged / 2 + T_flange / 2, 2));
              var W_flange_segment = 2 * I_flange_segment / (weld_leg_flanged + T_flange);
              var A_flange_segment = 2 * S_cyl * weld_leg_flanged;

              var sigma_bending = M_flange / W_flange_segment;
              this.set("sigma_bending", sigma_bending);

              var tau_shear = F_per_segment / A_flange_segment;
              this.set("tau_shear", tau_shear);

              var sigma_eqi = Math.sqrt(Math.pow(sigma_bending, 2) + 3 * Math.pow(tau_shear, 2));
              sigma_eqi = roundToDecimal(sigma_eqi, 3);
              this.set("sigma_eqi", sigma_eqi);

              if (sigma_eqi > sigma_flange)
              {
                this.error([], "The equivalent stress in the flange from pushing force is higher than \
                  the allowable stress: " + sigma_eqi + " MPa > " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }
              else
              {
                this.warn([], "The equivalent stress in the flange from pushing force is lower than or \
                  equal to the allowable stress: " + sigma_eqi + " MPa =< " + sigma_flange + " MPa. \
                  Please refer to $ref.");
              }

              var F_per_segment_2 = F_pull / n_bolt_flange;
              var M_flange_2 = F_per_segment_2 * L_moment;

              var sigma_bending_2 = M_flange_2 / (W_flange_segment);
              this.set("sigma_bending_2", sigma_bending_2);

              var tau_shear_2 = F_per_segment_2 / (A_flange_segment * 2);
              this.set("tau_shear_2", tau_shear_2);

              var sigma_eqi_2 = Math.sqrt(Math.pow(sigma_bending_2, 2) + 3 * Math.pow(tau_shear_2, 2));
              sigma_eqi_2 = roundToDecimal(sigma_eqi_2, 3);
              this.set("sigma_eqi_2", sigma_eqi_2);

              if (sigma_eqi_2 > sigma_flange)
              {
                this.error([], "The equivalent stress in the flange from pulling force is higher than the \
                  allowable stress: " + sigma_eqi_2 + " MPa > " + sigma_flange + " MPa. Please refer to $ref.");
              }
              else
              {
                this.warn([], "The equivalent stress in the flange from pulling force is lower than or equal the \
                  allowable stress: " + sigma_eqi_2 + " MPa =< " + sigma_flange + " MPa. Please refer to $ref.");
              }
              
              var sigma_eqi_fatigue_push = Math.sqrt(Math.pow(sigma_bending, 2) + 0.2 * Math.pow(tau_shear, 2));
              this.set("sigma_eqi_fatigue_push", sigma_eqi_fatigue_push);

              var sigma_eqi_fatigue_pull = Math.sqrt(Math.pow(sigma_bending_2, 2) + 0.2 * Math.pow(tau_shear_2, 2));
              this.set("sigma_eqi_fatigue_pull", sigma_eqi_fatigue_pull);

              var delta_stress = sigma_eqi_fatigue_push + sigma_eqi_fatigue_pull;
              this.set("delta_stress", delta_stress);

              var m_2 = 3;
              var log_a = 10.97;
              var log_N = log_a - m_2 * Math.log10(delta_stress);
              var n_fatigue = Math.pow(10, log_N);
              n_fatigue = Math.round(n_fatigue);

              if (n_fatigue < n_cycles_manufacturer_flanged) 
              {
                this.error([], "The flange's partial penetration weld weld has a predicted number of cycles to \
                  failure of: " + n_fatigue + " < " + n_cycles_manufacturer_flanged + ". Please refer to \
                  DNVGL-RP-C203.");
              }
              else
              {
                this.warn([], "The flange's partial penetration weld has a predicted number of cycles to \
                  failure of: " + n_fatigue + " >= " + n_cycles_manufacturer_flanged + ". Please refer to \
                  DNVGL-RP-C203.");
              }
            }
          }
        } // End of partial penetration flange weld rule
      ] // End of rules block
    } // End of app block
  }; // End of return
} // End of function
); // End of definition

