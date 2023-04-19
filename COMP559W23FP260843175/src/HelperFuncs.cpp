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
float clamp(float value, float min, float max) {
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
float max(float value1, float value2) {
	return (value1 > value2) ? value1 : value2;
}

/*
* Returns the lesser of both values
*/
float min(float value1, float value2) {
	return (value1 < value2) ? value1 : value2;
}
