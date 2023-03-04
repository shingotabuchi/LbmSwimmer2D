using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Runtime.InteropServices;

public class LbmSwimmer2DCompute : MonoBehaviour
{
    public Image plotImage;
    Texture2D plotTexture;
    Color[] plotPixels;
    ColorHeatMap colorHeatMap = new ColorHeatMap();
    public bool fixHeatMapMax;
    public bool fixHeatMapMin;
    public float fixedMaxSpeed = 0.6f;
    public float fixedMaxRho = 1.01f;
    public float fixedMinSpeed = 0f;
    public float fixedMinRho = 0f;
    public int DIM_X = 100;
    public int DIM_Y = 100;
    public int particleCount = 1;

    public float maxRho,minRho;
    public float maxSpeed,minSpeed;

    float[] cx = new float[9]{0, 1,    0,   -1,    0,     1,    -1,    -1,     1};
    float[] cy = new float[9]{0, 0,    1,    0,   -1,     1,     1,    -1,    -1};
    float[] w = new float[9]{4f/9f,1f/9f,1f/9f,1f/9f,1f/9f,1f/36f,1f/36f,1f/36f,1f/36f};
    float[] rho, u, v, speed, fx ,fy;
    float[] f, f0, ftmp;
    float tmp,u2,nu,tmp1, tmp2, tmp3,dt;
    public float particleDensity = 1.25f;
    public float particleRadius;
    public float epsw = 100.0f, zeta = 1.0f,tau = 0.53f;

    public int loopCount = 1;
    public float squirmerBeta = 0f;
    public float squirmerSpeedConstant = 0.001f;
    float[] gravity = new float[2];
    RoundParticle[] roundParticles;

    public ComputeShader compute;
    int immersedBoundary;
    ComputeBuffer roundParticlesSmallDataBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterPosBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterVelBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterFluidVelBuffer;
    ComputeBuffer roundParticleRoundParticleForceOnPerimeterBuffer;
    // Start is called before the first frame update
    void Start()
    {
        roundParticlesSmallDataBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticleSmallData)));
        roundParticleRoundParticlePerimeterPosBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterPos)));
        roundParticleRoundParticlePerimeterVelBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterVel)));
        roundParticleRoundParticlePerimeterFluidVelBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterFluidVel)));
        roundParticleRoundParticleForceOnPerimeterBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticleForceOnPerimeter)));

        immersedBoundary = compute.FindKernel("ImmersedBoundary");
        compute.SetBuffer(immersedBoundary,"roundParticlesSmallDataBuffer",roundParticlesSmallDataBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterPosBuffer",roundParticleRoundParticlePerimeterPosBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterVelBuffer",roundParticleRoundParticlePerimeterVelBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterFluidVelBuffer",roundParticleRoundParticlePerimeterFluidVelBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticleForceOnPerimeterBuffer",roundParticleRoundParticleForceOnPerimeterBuffer);

        // particleRadius = 0.125f * (float)DIM_X /2f;

        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotPixels = plotTexture.GetPixels();
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080 * (float)DIM_X/(float)DIM_Y,1080);
        
        speed = new float[DIM_X*DIM_Y];
        u = new float[DIM_X*DIM_Y];
        v = new float[DIM_X*DIM_Y];
        fx = new float[DIM_X*DIM_Y];
        fy = new float[DIM_X*DIM_Y];
        rho = new float[DIM_X*DIM_Y];
        f = new float[9*DIM_X*DIM_Y];
        f0 = new float[9*DIM_X*DIM_Y];
        ftmp = new float[9*DIM_X*DIM_Y];

        nu = (tau - 0.5f)/3.0f;

        dt = nu/((float)DIM_X/2.0f)/((float)DIM_X/2.0f)/0.1f;
        gravity[1] = -981.0f/2.0f*(float)(DIM_X)*dt*dt;
        gravity[0] = 0.0f;
        roundParticles = new RoundParticle[particleCount];

        for (int n = 0; n < particleCount; n++)
        {
            roundParticles[n] = new RoundParticle(particleDensity,particleRadius,new float[2]{DIM_X/2,DIM_Y/2});
        }
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u[i + DIM_X*j] = 0f; v[i + DIM_X*j] = 0.0f;
                fx[i + DIM_X*j] = 0.0f; fy[i + DIM_X*j] = 0.0f;
                rho[i + DIM_X*j] = 1.0f;
                u2 = u[i + DIM_X*j]*u[i + DIM_X*j] + v[i + DIM_X*j]*v[i + DIM_X*j];  
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + DIM_X*j] + cy[k]*v[i + DIM_X*j];  
                    f0[k + 9*(i + DIM_X*j)] = w[k]*rho[i + DIM_X*j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k + 9*(i + DIM_X*j)] = f0[k + 9*(i + DIM_X*j)];
                }
                speed[i + DIM_X*j] = Mathf.Sqrt(u2);
                maxSpeed = Mathf.Max(maxSpeed,speed[i + DIM_X*j]);
                minSpeed = Mathf.Min(minSpeed,speed[i + DIM_X*j]);
                maxRho = Mathf.Max(maxRho,rho[i + DIM_X*j]);
                maxRho = Mathf.Min(maxRho,rho[i + DIM_X*j]);
            } 
        }
    }

    void LBMStep()
    {
        Collision();
        Streaming();
        BouncebackBoundaries();
        UpdateSpeedAndDensity();

        // roundParticlesBuffer.SetData(roundParticles);
        // roundParticlesSmallDataBuffer.SetData(new RoundParticleSmallData[1]{roundParticles[0].smallData});
        // roundParticleRoundParticlePerimeterPosBuffer.SetData(new RoundParticlePerimeterPos[1]{roundParticles[0].perimeterPos});
        // roundParticleRoundParticlePerimeterVelBuffer.SetData(new RoundParticlePerimeterVel[1]{roundParticles[0].perimeterVel});
        // roundParticleRoundParticlePerimeterFluidVelBuffer.SetData(new RoundParticlePerimeterFluidVel[1]{roundParticles[0].perimeterFluidVel});
        // roundParticleRoundParticleForceOnPerimeterBuffer.SetData(new RoundParticleForceOnPerimeter[1]{roundParticles[0].forceOnPerimeter});
        ImmersedBoundary();
    }

    // Update is called once per frame
    void Update()
    {
        for (int i = 0; i < loopCount; i++)
        {
            LBMStep();
        }
        UpdatePlot();
    }

    void UpdatePlot()
    {
        if(maxSpeed == 0f)maxSpeed = 1f;
        if(fixHeatMapMax)
        {
            maxSpeed = fixedMaxSpeed;
        }
        if(fixHeatMapMin)
        {
            minSpeed = fixedMinSpeed;
        }
        for (int i = 0; i < plotPixels.Length; i++)
        {
            plotPixels[i] = colorHeatMap.GetColorForValue(speed[i%DIM_X + (i/DIM_X)*DIM_X]-minSpeed,maxSpeed-minSpeed);
        }
        // for (int n = 0; n < particleCount; n++)
        // {
        //     roundParticles[n].PlotParticleFill(ref plotPixels,DIM_X);
        // }
        plotTexture.SetPixels(plotPixels);
        plotTexture.Apply();
    }

    void Collision()
    {
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                u2 = u[i + DIM_X*j]*u[i + DIM_X*j] + v[i + DIM_X*j]*v[i + DIM_X*j];  
                for (int k = 0; k < 9; k++)
                {
                    tmp = cx[k]*u[i + DIM_X*j] + cy[k]*v[i + DIM_X*j];  
                    f0[k + 9*(i + DIM_X*j)] = w[k]*rho[i + DIM_X*j]*(1.0f +3.0f*tmp +9.0f/2.0f*tmp*tmp -3.0f/2.0f*u2);
                    f[k + 9*(i + DIM_X*j)] = f[k + 9*(i + DIM_X*j)] - (f[k + 9*(i + DIM_X*j)] - f0[k + 9*(i + DIM_X*j)])/tau + 3f*w[k]*(fx[i + DIM_X*j]*cx[k] + fy[i + DIM_X*j]*cy[k]);
                }
                // reset force;
                fx[i + DIM_X*j] = 0.0f;
                fy[i + DIM_X*j] = 0.0f;
            }   
        }
    }

    void Streaming()
    {
        ftmp = (float[])(f.Clone());

        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            { 
                for(int k = 0; k < 9; k++)
                {
                    // periodic boundary
                    int im = (i + (int)cx[k] + DIM_X)%DIM_X; 
                    int jm = (j + (int)cy[k] + DIM_Y)%DIM_Y;
                    // int im = i + (int)cx[k]; 
                    // int jm = j + (int)cy[k];
                    if((jm!=DIM_Y&&jm!=-1) && (im!=DIM_X&&im!=-1))
                    {
                        f[k + 9*(im + DIM_X*jm)] = ftmp[k + 9*(i + DIM_X*j)];
                    }
                } 
            }
        }
    }

    void BouncebackBoundaries()
    {
        // for (int i = 0; i < DIM_X; i++)
        // {
        //     f[4,i,DIM_Y-1] = f[2,i,DIM_Y-1];
        //     f[7,i,DIM_Y-1] = f[5,i,DIM_Y-1];
        //     f[8,i,DIM_Y-1] = f[6,i,DIM_Y-1]; 
        //     f[2,i,0] = f[4,i,0]; 
        //     f[5,i,0] = f[7,i,0]; 
        //     f[6,i,0] = f[8,i,0]; 
        // }
        // for (int j = 0; j < DIM_Y; j++)
        // {
        //     f[3,DIM_X-1,j] = f[1,DIM_X-1,j];
        //     f[6,DIM_X-1,j] = f[8,DIM_X-1,j];
        //     f[7,DIM_X-1,j] = f[5,DIM_X-1,j];
        //     f[1,0,j] = f[3,0,j]; 
        //     f[5,0,j] = f[7,0,j]; 
        //     f[8,0,j] = f[6,0,j]; 
        // }
    }

    void UpdateSpeedAndDensity()
    {
        maxSpeed = 0f;
        minSpeed = Mathf.Infinity;
        maxRho = 0f;
        minRho = Mathf.Infinity;
        for(int i = 0; i < DIM_X; i++)
        { 
            for(int j = 0; j < DIM_Y; j++)
            {
                rho[i + DIM_X*j] = f[0 + 9*(i + DIM_X*j)]; 
                u[i + DIM_X*j] = 0; v[i + DIM_X*j] = 0;
                for(int k = 1; k <= 8; k++)
                {
                    rho[i + DIM_X*j] = rho[i + DIM_X*j] + f[k + 9*(i + DIM_X*j)];
                    u[i + DIM_X*j] =   u[i + DIM_X*j] + f[k + 9*(i + DIM_X*j)]*cx[k];
                    v[i + DIM_X*j] =   v[i + DIM_X*j] + f[k + 9*(i + DIM_X*j)]*cy[k];
                } 
                u[i + DIM_X*j] = u[i + DIM_X*j]/rho[i + DIM_X*j];
                v[i + DIM_X*j] = v[i + DIM_X*j]/rho[i + DIM_X*j];
                speed[i + DIM_X*j] = Mathf.Sqrt(u[i + DIM_X*j]*u[i + DIM_X*j] + v[i + DIM_X*j]*v[i + DIM_X*j]);
                maxSpeed = Mathf.Max(maxSpeed,speed[i + DIM_X*j]);
                maxRho = Mathf.Max(maxRho,rho[i + DIM_X*j]);
                minSpeed = Mathf.Min(minSpeed,speed[i + DIM_X*j]);
                minRho = Mathf.Min(minRho,rho[i + DIM_X*j]);
            } 
        }
    }

    void ImmersedBoundary()
    {
        for(int n = 0; n < particleCount; n++) 
        { 
            roundParticles[n].smallData.forceFromCollisions[0] = 0f;
            roundParticles[n].smallData.forceFromCollisions[1] = 0f;
            tmp1 = Mathf.Abs(roundParticles[n].smallData.pos[1] + roundParticles[n].smallData.radius); 
            if(tmp1 < 2.0f*roundParticles[n].smallData.radius + zeta){
                roundParticles[n].smallData.forceFromCollisions[1] = (roundParticles[n].smallData.pos[1] + roundParticles[n].smallData.radius)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(DIM_Y-1-roundParticles[n].smallData.pos[1] + roundParticles[n].smallData.radius); 
            if(tmp1 < 2.0f*roundParticles[n].smallData.radius + zeta){
                roundParticles[n].smallData.forceFromCollisions[1] = -(DIM_Y-1-roundParticles[n].smallData.pos[1] + roundParticles[n].smallData.radius)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(roundParticles[n].smallData.pos[0] + roundParticles[n].smallData.radius); 
            if(tmp1 < 2.0f*roundParticles[n].smallData.radius + zeta){
                roundParticles[n].smallData.forceFromCollisions[0] = (roundParticles[n].smallData.pos[0] + roundParticles[n].smallData.radius)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)/epsw;
            }
            tmp1 = Mathf.Abs(DIM_X-1-roundParticles[n].smallData.pos[0] + roundParticles[n].smallData.radius); 
            if(tmp1 < 2.0f*roundParticles[n].smallData.radius + zeta){
                roundParticles[n].smallData.forceFromCollisions[0] = -(DIM_X-1-roundParticles[n].smallData.pos[0] + roundParticles[n].smallData.radius)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)/epsw;
            }

            for (int k = 0; k < particleCount; k++)
            {
                if(k==n) continue;
                for (int i = 0; i < 2; i++)
                {
                    tmp1 = roundParticles[n].ParticleDistance(roundParticles[k]);
                    if(tmp1 < 2.0f*roundParticles[n].smallData.radius + zeta){
                        roundParticles[n].smallData.forceFromCollisions[i] += (roundParticles[n].smallData.pos[i] - roundParticles[k].smallData.pos[i])*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)*(2.0f*roundParticles[n].smallData.radius - tmp1 + zeta)/epsw;
                    }
                }
            }

            roundParticles[n].smallData.forceFromFluid[0] = 0f;
            roundParticles[n].smallData.forceFromFluid[1] = 0f;
            roundParticles[n].smallData.torque = 0f;
            for(int m = 0; m < roundParticles[n].smallData.perimeterPointCount ; m++) 
            {
                roundParticles[n].perimeterFluidVel.perimeterFluidVel[0 +2*m] = 0f;
                roundParticles[n].perimeterFluidVel.perimeterFluidVel[1 +2*m] = 0f;
                // 固体表面の速度を計算
                for(int i = (int)roundParticles[n].perimeterPos.perimeterPos[0 +2*m] - 3; i < (int)roundParticles[n].perimeterPos.perimeterPos[0 +2*m] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos.perimeterPos[1 +2*m] - 3; j < (int)roundParticles[n].perimeterPos.perimeterPos[1 +2*m] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos.perimeterPos[0 +2*m] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos.perimeterPos[1 +2*m] - (float)j);
                        if(tmp1 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if(tmp2 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            roundParticles[n].perimeterFluidVel.perimeterFluidVel[0 +2*m] += u[i + DIM_X*j]*tmp3;
                            roundParticles[n].perimeterFluidVel.perimeterFluidVel[1 +2*m] += v[i + DIM_X*j]*tmp3;
                        }
                    } 
                }
                float boundaryPointTheta = roundParticles[n].smallData.theta + m * 2f * Mathf.PI / roundParticles[n].smallData.perimeterPointCount;
                float sin = Mathf.Sin(boundaryPointTheta);
                float cos = Mathf.Cos(boundaryPointTheta);
                float surfaceVelocityNorm = squirmerSpeedConstant * ( sin + 2f * squirmerBeta * sin * cos );
                Vector2 surfaceVelocity = new Vector2(cos,-sin) * surfaceVelocityNorm;
                roundParticles[n].perimeterFluidVel.perimeterFluidVel[0 +2*m] += surfaceVelocity[0];
                roundParticles[n].perimeterFluidVel.perimeterFluidVel[1 +2*m] += surfaceVelocity[1];
                roundParticles[n].forceOnPerimeter.forceOnPerimeter[0 +2*m] = roundParticles[n].perimeterVel.perimeterVel[0 +2*m] - roundParticles[n].perimeterFluidVel.perimeterFluidVel[0 +2*m];
                roundParticles[n].forceOnPerimeter.forceOnPerimeter[1 +2*m] = roundParticles[n].perimeterVel.perimeterVel[1 +2*m] - roundParticles[n].perimeterFluidVel.perimeterFluidVel[1 +2*m];

                // 固体が外部に与える力を計算
                for(int i = (int)roundParticles[n].perimeterPos.perimeterPos[0 +2*m] - 3; i < (int)roundParticles[n].perimeterPos.perimeterPos[0 +2*m] + 3; i++)
                {
                    for(int j = (int)roundParticles[n].perimeterPos.perimeterPos[1 +2*m] - 3; j < (int)roundParticles[n].perimeterPos.perimeterPos[1 +2*m] + 3; j++)
                    {
                        tmp1 = Mathf.Abs(roundParticles[n].perimeterPos.perimeterPos[0 +2*m] - (float)i);
                        tmp2 = Mathf.Abs(roundParticles[n].perimeterPos.perimeterPos[1 +2*m] - (float)j);
                        if(tmp1 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp1/2.0f))/4.0f;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if(tmp2 <= 2.0f)
                        {
                            tmp3 = (1.0f + Mathf.Cos(Mathf.PI*tmp2/2.0f))/4.0f*tmp3;
                        } 
                        else 
                        {
                            tmp3 = 0.0f;
                        }
                        if((j<DIM_Y&&j>=0) && (i<DIM_X&&i>=0))
                        {
                            fx[i + DIM_X*j] += roundParticles[n].forceOnPerimeter.forceOnPerimeter[0 +2*m] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].smallData.radius/(float)roundParticles[n].smallData.perimeterPointCount;
                            fy[i + DIM_X*j] += roundParticles[n].forceOnPerimeter.forceOnPerimeter[1 +2*m] * tmp3 * 2.0f*Mathf.PI*roundParticles[n].smallData.radius/(float)roundParticles[n].smallData.perimeterPointCount;
                        }
                    } 
                }
                roundParticles[n].smallData.forceFromFluid[0] += roundParticles[n].forceOnPerimeter.forceOnPerimeter[0 +2*m];
                roundParticles[n].smallData.forceFromFluid[1] += roundParticles[n].forceOnPerimeter.forceOnPerimeter[1 +2*m];
                roundParticles[n].smallData.torque += roundParticles[n].forceOnPerimeter.forceOnPerimeter[1 +2*m] * (roundParticles[n].perimeterPos.perimeterPos[0 +2*m] - roundParticles[n].smallData.pos[0]) 
                                        - roundParticles[n].forceOnPerimeter.forceOnPerimeter[0 +2*m] * (roundParticles[n].perimeterPos.perimeterPos[1 +2*m] - roundParticles[n].smallData.pos[1]);
            } 

            roundParticles[n].smallData.forceFromFluid[0] *= -2f*Mathf.PI*roundParticles[n].smallData.radius/(float)roundParticles[n].smallData.perimeterPointCount;  
            roundParticles[n].smallData.forceFromFluid[1] *= -2f*Mathf.PI*roundParticles[n].smallData.radius/(float)roundParticles[n].smallData.perimeterPointCount;  
            roundParticles[n].smallData.torque *= -2f*Mathf.PI*roundParticles[n].smallData.radius/(float)roundParticles[n].smallData.perimeterPointCount;  

            roundParticles[n].UpdatePosVel();
            roundParticles[n].UpdateOmegaTheta();
            roundParticles[n].UpdatePerimeter();
        }
    }
}
