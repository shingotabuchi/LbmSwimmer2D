using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class LbmSwimmer2DProt : MonoBehaviour
{
    // public Image plotImage;
    // Texture2D plotTexture;
    // Color[] plotPixels;
    // ColorHeatMap colorHeatMap = new ColorHeatMap();
    // public bool fixHeatMapMax;
    // public bool fixHeatMapMin;
    // public float fixedMaxSpeed = 0.6f;
    // public float fixedMaxRho = 1.01f;
    // public float fixedMinSpeed = 0f;
    // public float fixedMinRho = 0f;
    // public int DIM_X = 100;
    // public int DIM_Y = 400;
    // public int particleCount = 2;

    // public float maxRho,minRho;
    // public float maxSpeed,minSpeed;

    // float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    // float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    // float[] w = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    // float[,] rho, u, v, speed, fx ,fy;
    // float[,,] f, f0, ftmp;
    // float tmp,u2,nu,tmp1, tmp2, tmp3,dt;
    // public float particleDensity = 1.25f;
    // float particleRadius;
    // public float epsw = 100.0f, zeta = 1.0f,tau = 0.53f;

    // public int loopCount = 1;
    // public float gRate = 6f;
    // float[] gravity = new float[2];
    // RoundParticle[] roundParticles;

    // // Start is called before the first frame update
    // void Start()
    // {
    //     particleRadius = 0.125f * (float)DIM_X /2f;

    //     plotTexture = new Texture2D(DIM_X,DIM_Y);
    //     plotTexture.filterMode = FilterMode.Point;
    //     plotPixels = plotTexture.GetPixels();
    //     plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
    //     ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080 * (float)DIM_X/(float)DIM_Y,1080);
        
    //     speed = new float[DIM_X,DIM_Y];
    //     u = new float[DIM_X,DIM_Y];
    //     v = new float[DIM_X,DIM_Y];
    //     fx = new float[DIM_X,DIM_Y];
    //     fy = new float[DIM_X,DIM_Y];
    //     rho = new float[DIM_X,DIM_Y];
    //     f = new float[9,DIM_X,DIM_Y];
    //     f0 = new float[9,DIM_X,DIM_Y];
    //     ftmp = new float[9,DIM_X,DIM_Y];

    //     nu = (tau - 0.5f)/3.0f;

    //     // dt = nu/((float)(DIM_X-1)/2.0f)/((float)(DIM_X-1)/2.0f)/0.1f;
    //     dt = nu/((float)DIM_X/2.0f)/((float)DIM_X/2.0f)/0.1f;
    //     // gravity[1] = -gRate*981.0f*(float)(DIM_X-1)*dt*dt;
    //     // gravity[1] = -gRate*9.81e-5f;
    //     gravity[1] = -981.0f/2.0f*(float)(DIM_X)*dt*dt;
    //     gravity[0] = 0.0f;

    //     roundParticles = new RoundParticle[particleCount];
    //     for (int n = 0; n < particleCount; n++)
    //     {
    //         // roundParticles[n] = new RoundParticle(particleDensity,particleRadius,new float[2]{20f + (n%7)*10f,300f + (n/7)*10f});
    //         roundParticles[n] = new RoundParticle(particleDensity,particleRadius,new float[2]{(float)DIM_X / 2f, (float)DIM_X * 2f});
    //     }
    //     maxSpeed = 0f;
    //     minSpeed = Mathf.Infinity;
    //     maxRho = 0f;
    //     minRho = Mathf.Infinity;
    //     for(int i = 0; i < DIM_X; i++)
    //     { 
    //         for(int j = 0; j < DIM_Y; j++)
    //         {
    //             u[i,j] = 0f; v[i,j] = 0.0f;
    //             fx[i,j] = 0.0f; fy[i,j] = 0.0f;
    //             rho[i,j] = 1.0f;
    //             u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];  
    //             for (int k = 0; k < 9; k++)
    //             {
    //                 tmp = cx[k]*u[i,j] + cy[k]*v[i,j];  
    //                 f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
    //                 f[k,i,j] = f0[k,i,j];
    //             }
    //             speed[i,j] = Mathf.Sqrt(u2);
    //             maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
    //             minSpeed = Mathf.Min(minSpeed,speed[i,j]);
    //             maxRho = Mathf.Max(maxRho,rho[i,j]);
    //             maxRho = Mathf.Min(maxRho,rho[i,j]);
    //         } 
    //     }
    // }

    // void LBMStep()
    // {
    //     Collision();
    //     Streaming();
    //     BouncebackBoundaries();
    //     UpdateSpeedAndDensity();
    //     ImmersedBoundary();
    // }

    // // Update is called once per frame
    // void Update()
    // {
    //     for (int i = 0; i < loopCount; i++)
    //     {
    //         LBMStep();
    //     }
    //     UpdatePlot();
    // }

    // void UpdatePlot()
    // {
    //     if(maxSpeed == 0f)maxSpeed = 1f;
    //     if(maxRho == 0f)maxRho = 1f;
    //     if(fixHeatMapMax)
    //     {
    //         maxSpeed = fixedMaxSpeed;
    //         maxRho = fixedMaxRho;
    //     }
    //     if(fixHeatMapMin)
    //     {
    //         minSpeed = fixedMinSpeed;
    //         minSpeed = fixedMinRho;
    //     }
    //     for (int i = 0; i < plotPixels.Length; i++)
    //     {
    //         plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X,i/DIM_X]-minSpeed,maxSpeed-minSpeed);
    //     }
    //     for (int n = 0; n < particleCount; n++)
    //     {
    //         roundParticles[n].PlotParticleFill(ref plotPixels,DIM_X);
    //     }
    //     plotTexture.SetPixels(plotPixels);
    //     plotTexture.Apply();
    // }

    // void Collision()
    // {
    //     for(int i = 0; i < DIM_X; i++)
    //     { 
    //         for(int j = 0; j < DIM_Y; j++)
    //         {
    //             u2 = u[i,j]*u[i,j] + v[i,j]*v[i,j];  
    //             for (int k = 0; k < 9; k++)
    //             {
    //                 tmp = cx[k]*u[i,j] + cy[k]*v[i,j];  
    //                 f0[k,i,j] = w[k]*rho[i,j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
    //                 f[k,i,j] = f[k,i,j] - (f[k,i,j] - f0[k,i,j])/tau + 3f*w[k]*(fx[i,j]*cx[k] + fy[i,j]*cy[k]);
    //             }
    //             // reset force;
    //             fx[i,j] = 0.0f;
    //             fy[i,j] = 0.0f;
    //         }   
    //     }
    // }

    // void Streaming()
    // {
    //     ftmp = (float[,,])(f.Clone());

    //     for(int i = 0; i < DIM_X; i++)
    //     { 
    //         for(int j = 0; j < DIM_Y; j++)
    //         { 
    //             for(int k = 0; k < 9; k++)
    //             {
    //                 // // periodic boundary
    //                 // int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
    //                 // int jm = (j + (int)cy[k] + DIM_Y)%DIM_Y;
    //                 int im = i + (int)cx[k]; 
    //                 int jm = j + (int)cy[k];
    //                 if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
    //                 {
    //                     f[k,im,jm] = ftmp[k,i,j];
    //                 }
    //             } 
    //         }
    //     }
    // }

    // void BouncebackBoundaries()
    // {
    //     for (int i = 0; i < DIM_X; i++)
    //     {
    //         f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
    //         f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
    //         f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
    //         f[2,i,0] = f[4,i,0]; 
    //         f[5,i,0] = f[7,i,0]; 
    //         f[6,i,0] = f[8,i,0]; 
    //     }
    //     for (int j = 0; j < DIM_Y; j++)
    //     {
    //         f[3,DIM_X-1,j] = f[1,DIM_X-1,j];
    //         f[6,DIM_X-1,j] = f[8,DIM_X-1,j];
    //         f[7,DIM_X-1,j] = f[5,DIM_X-1,j];
    //         f[1,0,j] = f[3,0,j]; 
    //         f[5,0,j] = f[7,0,j]; 
    //         f[8,0,j] = f[6,0,j]; 
    //     }
    // }

    // void UpdateSpeedAndDensity()
    // {
    //     maxSpeed = 0f;
    //     minSpeed = Mathf.Infinity;
    //     maxRho = 0f;
    //     minRho = Mathf.Infinity;
    //     for(int i = 0; i < DIM_X; i++)
    //     { 
    //         for(int j = 0; j < DIM_Y; j++)
    //         {
    //             rho[i,j] = f[0,i,j]; 
    //             u[i,j] = 0; v[i,j] = 0;
    //             for(int k = 1; k <= 8; k++)
    //             {
    //                 rho[i,j] = rho[i,j] + f[k,i,j];
    //                 u[i,j] =   u[i,j] + f[k,i,j]*cx[k];
    //                 v[i,j] =   v[i,j] + f[k,i,j]*cy[k];
    //             } 
    //             u[i,j] = u[i,j]/rho[i,j];
    //             v[i,j] = v[i,j]/rho[i,j];
    //             speed[i,j] = Mathf.Sqrt(u[i,j]*u[i,j] + v[i,j]*v[i,j]);
    //             maxSpeed = Mathf.Max(maxSpeed,speed[i,j]);
    //             maxRho = Mathf.Max(maxRho,rho[i,j]);
    //             minSpeed = Mathf.Min(minSpeed,speed[i,j]);
    //             minRho = Mathf.Min(minRho,rho[i,j]);
    //         } 
    //     }
    // }

    // void ImmersedBoundary()
    // {
    //     for(int n = 0; n < particleCount; n++) 
    //     { 
    //         roundParticles[n].forceFromCollisions[0] = 0f;
    //         roundParticles[n].forceFromCollisions[1] = 0f;
    //         tmp1 = Mathf.Abs(roundParticles[n].pos[1] + roundParticles[n].radius); 
    //         if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
    //             roundParticles[n].forceFromCollisions[1] = (roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
    //         }
    //         tmp1 = Mathf.Abs(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius); 
    //         if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
    //             roundParticles[n].forceFromCollisions[1] = -(DIM_Y-1-roundParticles[n].pos[1] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
    //         }
    //         tmp1 = Mathf.Abs(roundParticles[n].pos[0] + roundParticles[n].radius); 
    //         if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
    //             roundParticles[n].forceFromCollisions[0] = (roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
    //         }
    //         tmp1 = Mathf.Abs(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius); 
    //         if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
    //             roundParticles[n].forceFromCollisions[0] = -(DIM_X-1-roundParticles[n].pos[0] + roundParticles[n].radius)*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
    //         }

    //         for (int k = 0; k < particleCount; k++)
    //         {
    //             if(k==n) continue;
    //             for (int i = 0; i < 2; i++)
    //             {
    //                 tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
    //                 if(tmp1 < 2.0f*roundParticles[n].radius + zeta){
    //                     roundParticles[n].forceFromCollisions[i] += (roundParticles[n].pos[i] - roundParticles[k].pos[i])*(2.0f*roundParticles[n].radius - tmp1 + zeta)*(2.0f*roundParticles[n].radius - tmp1 + zeta)/epsw;
    //                 }
    //             }
    //         }

    //         roundParticles[n].forceFromFluid[0] = 0f;
    //         roundParticles[n].forceFromFluid[1] = 0f;
    //         roundParticles[n].torque = 0f;
    //         for(int m = 0; m < roundParticles[n].perimeterPointCount ; m++) 
    //         {
    //             roundParticles[n].perimeterFluidVel[m,0] = 0f;
    //             roundParticles[n].perimeterFluidVel[m,1] = 0f;
    //             // 固体表面の速度を計算
    //             for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
    //             {
    //                 for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
    //                 {
    //                     tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
    //                     tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
    //                     if(tmp1 <= 2.0f)
    //                     {
    //                         tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
    //                     } 
    //                     else 
    //                     {
    //                         tmp3 = 0.0f;
    //                     }
    //                     if(tmp2 <= 2.0f)
    //                     {
    //                         tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
    //                     } 
    //                     else 
    //                     {
    //                         tmp3 = 0.0f;
    //                     }
    //                     if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
    //                     {
    //                         roundParticles[n].perimeterFluidVel[m,0] += u[i,j]*tmp3;
    //                         roundParticles[n].perimeterFluidVel[m,1] += v[i,j]*tmp3;
    //                     }
    //                 } 
    //             }
    //             roundParticles[n].forceOnPerimeter[m,0] = roundParticles[n].perimeterVel[m,0] - roundParticles[n].perimeterFluidVel[m,0];
    //             roundParticles[n].forceOnPerimeter[m,1] = roundParticles[n].perimeterVel[m,1] - roundParticles[n].perimeterFluidVel[m,1];

    //             // 固体が外部に与える力を計算
    //             for(int i = (int)roundParticles[n].perimeterPos[m,0] - 3; i < (int)roundParticles[n].perimeterPos[m,0] + 3; i++)
    //             {
    //                 for(int j = (int)roundParticles[n].perimeterPos[m,1] - 3; j < (int)roundParticles[n].perimeterPos[m,1] + 3; j++)
    //                 {
    //                     tmp1 = Mathf.Abs(roundParticles[n].perimeterPos[m,0] - (float)i);
    //                     tmp2 = Mathf.Abs(roundParticles[n].perimeterPos[m,1] - (float)j);
    //                     if(tmp1 <= 2.0f)
    //                     {
    //                         tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
    //                     } 
    //                     else 
    //                     {
    //                         tmp3 = 0.0f;
    //                     }
    //                     if(tmp2 <= 2.0f)
    //                     {
    //                         tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
    //                     } 
    //                     else 
    //                     {
    //                         tmp3 = 0.0f;
    //                     }
    //                     if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
    //                     {
    //                         fx[i,j] += roundParticles[n].forceOnPerimeter[m,0] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
    //                         fy[i,j] += roundParticles[n].forceOnPerimeter[m,1] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;
    //                     }
    //                 } 
    //             }
    //             roundParticles[n].forceFromFluid[0] += roundParticles[n].forceOnPerimeter[m,0];
    //             roundParticles[n].forceFromFluid[1] += roundParticles[n].forceOnPerimeter[m,1];
    //             roundParticles[n].torque += roundParticles[n].forceOnPerimeter[m,1] * (roundParticles[n].perimeterPos[m,0] - roundParticles[n].pos[0]) 
    //                                     - roundParticles[n].forceOnPerimeter[m,0] * (roundParticles[n].perimeterPos[m,1] - roundParticles[n].pos[1]);
    //         } 

    //         roundParticles[n].forceFromFluid[0] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
    //         roundParticles[n].forceFromFluid[1] *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  
    //         roundParticles[n].torque *= -2f*Mathf.PI*roundParticles[n].radius/(float)roundParticles[n].perimeterPointCount;  

    //         roundParticles[n].UpdatePosVel(gravity);
    //         roundParticles[n].UpdateOmegaTheta();
    //         roundParticles[n].UpdatePerimeter();
    //     }
    // }
}
