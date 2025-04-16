#version 330 core
out vec4 FragColor;

in vec3 Normal;
in vec3 FragPos;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 objectColor;

void main() {
    // 增强环境光
    float ambientStrength = 0.5;  // 从0.3提高到0.5
    vec3 ambient = ambientStrength * lightColor;
    
    // 增强漫反射
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor * 0.8;  // 增加漫反射强度
    
    // 调整高光
    float specularStrength = 0.4;  // 降低高光强度，使表面不那么"塑料感"
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;
    
    // 确保基础亮度
    float minBrightness = 0.3;  // 添加最小亮度保证
    vec3 result = max((ambient + diffuse + specular) * objectColor, minBrightness * objectColor);
    
    FragColor = vec4(result, 1.0);

    result = (ambient + diffuse + specular) * objectColor;
    
    // 计算重心坐标到三角形边的距离
    float d = min(min(gl_FragCoord.x, gl_FragCoord.y), gl_FragCoord.z);
    float wireframe = 1.0 - smoothstep(0.0, 1.0, d);
    
    // 混合网格线
    vec3 wireframeColor = vec3(0.2);  // 网格线颜色
    result = mix(result, wireframeColor, wireframe * 0.5);  // 0.5是网格线强度
    
    FragColor = vec4(result, 1.0);
}