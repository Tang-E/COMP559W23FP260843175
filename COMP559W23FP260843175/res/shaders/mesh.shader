#shader vertex
attribute vec2 attrPosition;
uniform vec2 domainSize;
uniform vec3 color;
uniform vec2 translation;
uniform float scale;
varying vec3 fragColor;
void main() {
	vec2 v = translation + attrPosition * scale;
	vec4 screenTransform =
		vec4(2.0 / domainSize.x, 2.0 / domainSize.y, -1.0, -1.0);
	gl_Position =
		vec4(v * screenTransform.xy + screenTransform.zw, 0.0, 1.0);
	fragColor = color;
}

#shader fragment
precision mediump float;
varying vec3 fragColor;
void main() {
	gl_FragColor = vec4(fragColor, 1.0);
}
