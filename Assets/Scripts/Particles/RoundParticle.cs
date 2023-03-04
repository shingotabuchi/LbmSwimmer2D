using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;
public struct RoundParticleSmallData
{
    public float density;
    public float omega;//角速度
    public float theta;//角
    public float prevOmega1;//前フレームの角速度
    public float prevOmega2;//前々フレームの角速度
    public float torque;
    public int perimeterPointCount;
    public float volume;
    public float mass;
    public float momentOfInertia;
    public float radius;
    public int bottomPointIndex;

    public Vector2 pos;
    public Vector2 vel;
    public Vector2 prevVel1;//前フレームの速度
    public Vector2 prevVel2;//前々フレームの速度
    public Vector2 forceFromCollisions;
    public Vector2 forceFromFluid;
}

// [StructLayout(LayoutKind.Sequential)]
// public struct RoundParticlePerimeterPos
// {
//     [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2*252)]
//     public float[] perimeterPos;
// };
// [StructLayout(LayoutKind.Sequential)]
// public struct RoundParticlePerimeterVel
// {
//     [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2*252)]
//     public float[] perimeterVel;
// };
// [StructLayout(LayoutKind.Sequential)]
// public struct RoundParticlePerimeterFluidVel
// {
//     [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2*252)]
//     public float[] perimeterFluidVel;
// };
// [StructLayout(LayoutKind.Sequential)]
// public struct RoundParticleForceOnPerimeter
// {
//     [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2*252)]
//     public float[] forceOnPerimeter;
// };

// [StructLayout(LayoutKind.Sequential)]
// public struct RoundParticle
// {
//     public RoundParticleSmallData smallData;
//     public RoundParticlePerimeterPos perimeterPos;
//     public RoundParticlePerimeterVel perimeterVel;
//     public RoundParticlePerimeterFluidVel perimeterFluidVel;
//     public RoundParticleForceOnPerimeter forceOnPerimeter;
//     // public void PlotParticlePerimeter(ref Color[] pixels, int dim_x, Color? color = null)
//     // {
//     //     color ??= Color.white;
//     //     for(int i = 0; i < smallData.perimeterPointCount ; i++) 
//     //     {
//     //         if(
//     //             (int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x < pixels.Length
//     //             &&  
//     //             (int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x >= 0
//     //         )
//     //         pixels[(int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x] = (Color)color;
//     //     }
//     // }
//     // public void PlotParticleTrace(ref Color[] pixels, int dim_x, Color? color = null)
//     // {
//     //     color ??= Color.white;
//     //     for(int i = 0; i < smallData.perimeterPointCount ; i++) 
//     //     {
//     //         if(
//     //             (int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x < pixels.Length
//     //             &&  
//     //             (int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x >= 0
//     //         )
//     //         pixels[(int)perimeterPos.perimeterPos[0 + 2*i] + (int)perimeterPos.perimeterPos[1 + 2*i] * dim_x] = (Color)color;
//     //     }
//     // }
//     // public void UpdatePosVel(float[]? gravity = null)
//     // {
//     //     gravity ??= new float[2]{0f,0f};
//     //     for (int i = 0; i < 2; i++)
//     //     {
//     //         smallData.vel[i] = (1f + 1f/smallData.density) * smallData.prevVel1[i]
//     //                                 - 1f/smallData.density * smallData.prevVel2[i]
//     //                                 + (smallData.forceFromFluid[i] + smallData.forceFromCollisions[i])/smallData.mass
//     //                                 + (1f - 1f/smallData.density) * (float)gravity[i];
//     //         smallData.pos[i] += (smallData.vel[i] + smallData.prevVel1[i])/2f;
//     //         smallData.prevVel2[i] = smallData.prevVel1[i];
//     //         smallData.prevVel1[i] = smallData.vel[i];
//     //     }
//     // }
//     // public void UpdateOmegaTheta()
//     // {
//     //     smallData.omega = (1f + 1f/smallData.density) * smallData.prevOmega1 
//     //                             - 1f/smallData.density * smallData.prevOmega2
//     //                             + smallData.torque/smallData.momentOfInertia;
//     //     smallData.theta += (smallData.omega + smallData.prevOmega1)/2f;
//     //     smallData.prevOmega2 = smallData.prevOmega1;
//     //     smallData.prevOmega1 = smallData.omega;
//     // }
    
//     // public RoundParticle(float _density,float _radius,float[] _initPos, float _initTheta = 0f)
//     // {
//     //     smallData.density = _density;
//     //     smallData.radius = _radius;
//     //     smallData.volume = Mathf.PI*smallData.radius*smallData.radius;
//     //     smallData.mass = smallData.volume * smallData.density;
//     //     smallData.momentOfInertia = (Mathf.PI * smallData.radius*smallData.radius*smallData.radius*smallData.radius * smallData.density)/2f;
//     //     smallData.pos = _initPos;
//     //     smallData.vel = new float[2]{0f,0f};
//     //     smallData.prevVel1 = new float[2]{0f,0f};
//     //     smallData.prevVel2 = new float[2]{0f,0f};
//     //     smallData.omega = 0f;
//     //     smallData.theta = _initTheta;
//     //     smallData.prevOmega1 = 0f;
//     //     smallData.prevOmega2 = 0f;
//     //     // max 252
//     //     smallData.perimeterPointCount = (int)(2.0f * Mathf.PI * smallData.radius * 2.0f);
//     //     if(smallData.perimeterPointCount > 252) Debug.Log("bruh");
//     //     // perimeterPos.perimeterPos = new float[smallData.perimeterPointCount*2];
//     //     // perimeterVel.perimeterVel = new float[smallData.perimeterPointCount*2];
//     //     // perimeterFluidVel.perimeterFluidVel = new float[smallData.perimeterPointCount*2];
//     //     // forceOnPerimeter.forceOnPerimeter = new float[smallData.perimeterPointCount*2];
//     //     smallData.forceFromCollisions = new float[2]{0f,0f};
//     //     smallData.forceFromFluid = new float[2]{0f,0f};
//     //     smallData.torque = 0f;

//     //     for(int i = 0; i < smallData.perimeterPointCount; i++) 
//     //     {
//     //         float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)smallData.perimeterPointCount + smallData.theta;

//     //         perimeterPos.perimeterPos[0 + 2*i] = smallData.pos[0] + smallData.radius * Mathf.Cos(angle);
//     //         perimeterPos.perimeterPos[1 + 2*i] = smallData.pos[1] + smallData.radius * Mathf.Sin(angle);
//     //         perimeterVel.perimeterVel[0 + 2*i] = smallData.vel[0] - smallData.omega*(perimeterPos.perimeterPos[1 + 2*i] - smallData.pos[1]);
//     //         perimeterVel.perimeterVel[1 + 2*i] = smallData.vel[1] + smallData.omega*(perimeterPos.perimeterPos[0 + 2*i] - smallData.pos[0]);
//     //         forceOnPerimeter.forceOnPerimeter[0 + 2*i] = 0f;
//     //         forceOnPerimeter.forceOnPerimeter[1 + 2*i] = 0f;
//     //         perimeterFluidVel.perimeterFluidVel[0 + 2*i] = 0f;
//     //         perimeterFluidVel.perimeterFluidVel[1 + 2*i] = 0f;
//     //     } 
//     // }

//     // public void PlotParticleFill(ref Color[] pixels, int dim_x, Color? color = null)
//     // {
//     //     color ??= Color.white;
//     //     for (int i = -((int)smallData.radius); i <= (int)smallData.radius; i++)
//     //     {
//     //         for (int j = -((int)smallData.radius); j <= (int)smallData.radius; j++)
//     //         {
//     //             if(i*i + j*j <= smallData.radius*smallData.radius)
//     //             {
//     //                 if(
//     //                     i + (int)smallData.pos[0] + (j + (int)smallData.pos[1])* dim_x < pixels.Length
//     //                     &&  
//     //                     i + (int)smallData.pos[0] + (j + (int)smallData.pos[1])* dim_x >= 0
//     //                 )
//     //                 pixels[i + (int)smallData.pos[0] + (j + (int)smallData.pos[1])* dim_x] = (Color)color;
//     //             }
//     //         }
//     //     }
//     // }

//     // public void UpdatePerimeter()
//     // {
//     //     for(int i = 0; i < smallData.perimeterPointCount; i++) 
//     //     {
//     //         float angle = (3f*Mathf.PI)/2f + 2.0f*Mathf.PI*(float)i/(float)smallData.perimeterPointCount + smallData.theta;;
//     //         perimeterPos.perimeterPos[0 + 2*i] = smallData.pos[0] + smallData.radius * Mathf.Cos(angle);
//     //         perimeterPos.perimeterPos[1 + 2*i] = smallData.pos[1] + smallData.radius * Mathf.Sin(angle);
//     //         perimeterVel.perimeterVel[0 + 2*i] = smallData.vel[0] - smallData.omega*(perimeterPos.perimeterPos[1 + 2*i] - smallData.pos[1]);
//     //         perimeterVel.perimeterVel[1 + 2*i] = smallData.vel[1] + smallData.omega*(perimeterPos.perimeterPos[0 + 2*i] - smallData.pos[0]);
//     //     } 
//     // }

//     // public float ParticleDistance(RoundParticle particle)
//     // {
//     //     return Mathf.Sqrt( (smallData.pos[0]-particle.smallData.pos[0])*(smallData.pos[0]-particle.smallData.pos[0]) + (smallData.pos[1]-particle.smallData.pos[1])*(smallData.pos[1]-particle.smallData.pos[1]) );
//     // }

//     // public bool PointIsInParticle(int[] point)
//     // {
//     //     if(
//     //         (point[0] - smallData.pos[0])*(point[0] - smallData.pos[0]) + (point[1] - smallData.pos[1])*(point[1] - smallData.pos[1])
//     //         < smallData.radius*smallData.radius
//     //     )return true;
//     //     return false;
//     // }
// }
