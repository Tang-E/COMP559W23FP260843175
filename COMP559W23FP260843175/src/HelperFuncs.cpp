/*
* A bunch of handy helper functions that should kind of
* just exist but not be a particular method of any of the
* classes.
* 
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
*/

/*
* Clamps value between min and max.
*/
double clamp(double value, double min, double max) {
	if (value < min) {
		return min;
	}
	else if (value > max) {
		return max;
	}
	else {
		return value;
	}
}

/*
* Returns the greater of both values
*/
double max(double value1, double value2) {
	return (value1 > value2) ? value1 : value2;
}

/*
* Returns the lesser of both values
*/
double min(double value1, double value2) {
	return (value1 < value2) ? value1 : value2;
}

/*
* Provided the two pointers for two double arrays and
* their lengths (LENGTHS MUST MATCH), copies entries
* from the source array to the destination array.
*/
void copy(double* srcArr, double* destArr, int length) {
	for (int i = 0; i < length; i++) {
		destArr[i] = srcArr[i];
	}
}
