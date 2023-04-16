using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Runtime.InteropServices;
using UnityEngine.VFX;
public class ShanChen : MonoBehaviour
{
    public enum PlotMode
    {
        Phase,
        Density,
        Speed,
    }
    public PlotMode plotMode;
    public Image plotImage;
    Texture2D plotTexture;
    RenderTexture renderTexture;
    public int lx,ly,lz;
    public int plotXIndex;
    public int loopCount;
    public bool debugMode;
    public bool nextFrame;
    public ComputeShader compute;
    ComputeBuffer obst,uvw,psx,rho,ff;

    public float minRho,maxRho;
    public float rhow = 0.12f;
    public float tauc = 1.0f;
    public int eosIndex = 5;
    public float TT0,rho_h,rho_l;
    float[] TT0W = new float[]{0.975f, 0.95f, 0.925f, 0.9f, 0.875f, 0.85f, 0.825f, 0.8f, 0.775f, 0.75f, 0.7f, 0.65f};
	float[] RHW = new float[]{0.16f, 0.21f, 0.23f, 0.247f, 0.265f, 0.279f, 0.29f, 0.314f, 0.30f, 0.33f, 0.36f, 0.38f};
	float[] RLW = new float[]{0.08f, 0.067f, 0.05f, 0.0405f, 0.038f, 0.032f, 0.025f, 0.0245f, 0.02f, 0.015f, 0.009f, 0.006f};
    float[] rhoDebug,uvwDebug,psxDebug,ffDebug;
    int[] obstDebug;
    int init,plotDensity;
    int stream,calcPressureAndUvw,collision;

    private void OnDestroy() {
        obst.Dispose();
        uvw.Dispose();
        psx.Dispose();
        rho.Dispose();
        ff.Dispose();
    }
    int stepCount = 0;
    private void Start() 
    {
        plotTexture = new Texture2D(ly,lz);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,ly,lz),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(ly*1080/lz,1080);

        renderTexture = new RenderTexture(ly,lz,24);
        renderTexture.enableRandomWrite = true;

        obst = new ComputeBuffer(lx*ly*lz,sizeof(int));
        obstDebug = new int[lx*ly*lz];
        uvw = new ComputeBuffer(lx*ly*lz*3,sizeof(float));
        uvwDebug = new float[lx*ly*lz*3];
        rho = new ComputeBuffer(lx*ly*lz,sizeof(float));
        rhoDebug = new float[lx*ly*lz];
        psx = new ComputeBuffer(lx*ly*lz,sizeof(float));
        psxDebug = new float[lx*ly*lz];
        ff = new ComputeBuffer(lx*ly*lz*19*2,sizeof(float));
        ffDebug = new float[lx*ly*lz*19*2];

        compute.SetInt("lx",lx);
        compute.SetInt("ly",ly);
        compute.SetInt("lz",lz);
        OnValidate();

        init = compute.FindKernel("Init");
        compute.SetBuffer(init,"obst",obst);
        compute.SetBuffer(init,"uvw",uvw);
        compute.SetBuffer(init,"rho",rho);
        compute.SetBuffer(init,"ff",ff);

        plotDensity = compute.FindKernel("PlotDensity");
        compute.SetBuffer(plotDensity,"rho",rho);
        compute.SetTexture(plotDensity,"renderTexture",renderTexture);

        stream = compute.FindKernel("Stream");
        compute.SetBuffer(stream,"ff",ff);

        calcPressureAndUvw = compute.FindKernel("CalcPressureAndUvw");
        compute.SetBuffer(calcPressureAndUvw,"obst",obst);
        compute.SetBuffer(calcPressureAndUvw,"ff",ff);
        compute.SetBuffer(calcPressureAndUvw,"uvw",uvw);
        compute.SetBuffer(calcPressureAndUvw,"rho",rho);
        compute.SetBuffer(calcPressureAndUvw,"psx",psx);

        collision = compute.FindKernel("Collision");
        compute.SetBuffer(collision,"obst",obst);
        compute.SetBuffer(collision,"ff",ff);
        compute.SetBuffer(collision,"uvw",uvw);
        compute.SetBuffer(collision,"rho",rho);
        compute.SetBuffer(collision,"psx",psx);

        Init();
        // rho.GetData(rhoDebug);
        using (StreamWriter writer = new StreamWriter(stepCount.ToString() + ".txt"))
        {
            for (int k = 0; k < lz; k++)
            {
                for (int j = 0; j < ly; j++)
                {
                    for (int i = 0; i < lx; i++)
                    {
                        writer.WriteLine(i.ToString() + " " + j.ToString() + " " + k.ToString() + " " + rhoDebug[i + j*lx + k*lx*ly].ToString("0.00000000"));
                    }
                }
            }
        }
        PlotDensity();
    }
    
    private void FixedUpdate() {
        for (int kk = 0; kk < loopCount; kk++)
        {
            if(nextFrame || !debugMode)
            {
                Stream();
                CalcPressureAndUvw();
                Collision();
                nextFrame = false;
                // rho.GetData(rhoDebug);
                DebugStep();
                stepCount++;
                using (StreamWriter writer = new StreamWriter(stepCount.ToString() + ".txt"))
                {
                    for (int k = 0; k < lz; k++)
                    {
                        for (int j = 0; j < ly; j++)
                        {
                            for (int i = 0; i < lx; i++)
                            {
                                writer.WriteLine(i.ToString() + " " + j.ToString() + " " + k.ToString() + " " + rhoDebug[i + j*lx + k*lx*ly].ToString("0.00000000"));
                            }
                        }
                    }
                }
            }
        }
        PlotDensity();
    }

    void Stream()
    {
        compute.Dispatch(stream,(lx+7)/8,(ly+7)/8,(lz+7)/8);
    }
    void CalcPressureAndUvw()
    {
        compute.Dispatch(calcPressureAndUvw,(lx+7)/8,(ly+7)/8,(lz+7)/8);
    }
    void Collision()
    {
        compute.Dispatch(collision,(lx+7)/8,(ly+7)/8,(lz+7)/8);
    }

    void Init()
    {
        compute.Dispatch(init,(lx+7)/8,(ly+7)/8,(lz+7)/8);

        for (int k = 0; k < lz; k++)
        {
            for (int j = 0; j < ly; j++)
            {
                for (int i = 0; i < lx; i++)
                {
                    InitDebug(i,j,k);
                }
            }
        }
    }

    void PlotDensity()
    {
        compute.Dispatch(plotDensity,(ly+7)/8,(lz+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
    }

    private void OnValidate() {
        compute.SetFloat("minRho",minRho);
        compute.SetFloat("maxRho",maxRho);
        TT0 = TT0W[eosIndex];
        rho_h = RHW[eosIndex];
        rho_l = RLW[eosIndex];
        // rho_h = 1f;
        // rho_l = 1f;
        // rhow = 1f;
        compute.SetFloat("rho_h",rho_h);
        compute.SetFloat("rho_l",rho_l);
        compute.SetFloat("rhow",rhow);
        compute.SetFloat("tauc",tauc);
        compute.SetInt("plotXIndex",plotXIndex);
    }

    int BufferIndex(int x, int y, int z)
    {
        return x + y*lx + z*lx*ly;
    }

    int BufferIndexUvw(int bufferIndex, int uvwIndex)
    {
        return bufferIndex * 3 + uvwIndex;
    }

    int BufferIndexK(int bufferIndex, int k)
    {
        return bufferIndex * 19 + k;
    }

    int BufferIndexKTmp(int bufferIndex, int k)
    {
        return lx*ly*lz*19 + BufferIndexK(bufferIndex,k);
    }
    int UMod(int n, int m)
    {
        return (int)((uint)(n + m)%(uint)m);
    }

    float cc = 1.0f;
    float c_squ = 1.0f/3.0f;
    float R = 1.0f;
    float b = 4.0f;
    float a = 1.0f;
    float[] t_k = new float[]{
        1.0f / 3.0f,
        1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f
    };
    float[] xc = new float[]{
        0.0f,
        1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        1.0f, 1.0f, -1.0f, -1.0f, 1.0f, -1.0f,
        1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    float[] yc = new float[]{
        0.0f,
        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        1.0f, -1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 1.0f, -1.0f, -1.0f
    };
    float[] zc = new float[]{
        0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, -1.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f,
        -1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f
    };
    int[] ex = new int[]{
        0,
        1, -1, 0, 0, 0, 0,
        1, 1, -1, -1, 1, -1,
        1, -1, 0, 0, 0, 0
    };
    int[] ey = new int[]{
        0,
        0, 0, 1, -1, 0, 0,
        1, -1, 1, -1, 0, 0,
        0, 0, 1, 1, -1, -1
    };
    int[] ez = new int[]{
        0,
        0, 0, 0, 0, 1, -1,
        0, 0, 0, 0, 1, 1,
        -1, -1, 1, -1, 1, -1
    };
    int[] opp = new int[] {0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};
    float RR = 15.0f;

    void InitDebug(int x,int y,int z)
    {
        int k;
        int index = BufferIndex(x,y,z);
        obstDebug[index] = 0;
        if(z == 0) obstDebug[index] = 1;
        int indexU = BufferIndexUvw(index,0);
        int indexV = BufferIndexUvw(index,1);
        int indexW = BufferIndexUvw(index,2);
        uvwDebug[indexU] = 0.0f;
        uvwDebug[indexV] = 0.0f;
        uvwDebug[indexW] = 0.0f;
        float u = uvwDebug[indexU];
        float v = uvwDebug[indexV];
        float w = uvwDebug[indexW];
        rhoDebug[index] = rho_l;
        float xf = (float)(x-(int)((uint)lx/2));
        float yf = (float)(y-(int)((uint)ly/2));
        float zf = (float)z - RR/3.0f;
        if(xf*xf + yf*yf + zf*zf < RR*RR)
        {
            rhoDebug[index] = rho_h;
        }
        float u_squ = u*u + v*v + w*w;
        for(k = 0; k < 19; k++)
        {
            float u_n = xc[k] * u + yc[k] * v + zc[k] * w;
            float feq = 
            t_k[k]*rhoDebug[index] 
            *(
                cc*u_n/c_squ + (u_n*cc)*(u_n*cc)/(2.0f*c_squ*c_squ)
                - u_squ/(2.0f*c_squ))
            +t_k[k]*rhoDebug[index];
            
            ffDebug[BufferIndexK(index,k)] = feq;
        }
        // rho[index] = 0;
        // for(k = 0; k < 19; k++)
        // {
        //     int indexK = BufferIndexK(index,k);
        //     rho[index] += ffDebug[indexK];
        // }
    }

    void DebugStep()
    {
        for (int k = 0; k < lz; k++)
        {
            for (int j = 0; j < ly; j++)
            {
                for (int i = 0; i < lx; i++)
                {
                    StreamDebug(i,j,k);
                    CalcPressureAndUvwDebug(i,j,k);
                    CollisionDebug(i,j,k);
                }
            }
        }
    }

    void StreamDebug(int x, int y, int z)
    {
        int k;
        int index = BufferIndex(x,y,z);
        for(k = 0; k < 19; k++)
        {
            int xTmp = UMod(x + ex[k],lx);
            int yTmp = UMod(y + ey[k],ly);
            int zTmp = UMod(z + ez[k],lz);
            int indexTmp = BufferIndex(xTmp,yTmp,zTmp);
            ffDebug[BufferIndexKTmp(indexTmp,k)] = ffDebug[BufferIndexK(index,k)];
        }
    }

    void CalcPressureAndUvwDebug(int x, int y, int z)
    {
        float Tc = 0.3773f*a/(b*R);
        float TT= TT0 *Tc;
        int k;
        int index = BufferIndex(x,y,z);
        int indexU = BufferIndexUvw(index,0);
        int indexV = BufferIndexUvw(index,1);
        int indexW = BufferIndexUvw(index,2);
        
        uvwDebug[indexU] = 0.0f;
        uvwDebug[indexV] = 0.0f; 
        uvwDebug[indexW] = 0.0f; 
        rhoDebug[index] = 0.0f;

        for(k = 0; k < 19; k++)
        {
            int indexKTmp = BufferIndexKTmp(index,k);
            int indexK = BufferIndexK(index,k);
            ffDebug[indexK] = ffDebug[indexKTmp];
        }
        
        if(obstDebug[index] != 0)
        {
            // rhoDebug[index] = rhow;
            return;
        }

        for(k = 0; k < 19; k++)
        {
            int indexK = BufferIndexK(index,k);
            rhoDebug[index] += ffDebug[indexK];
            uvwDebug[indexU] += ffDebug[indexK] * xc[k];
            uvwDebug[indexV] += ffDebug[indexK] * yc[k];
            uvwDebug[indexW] += ffDebug[indexK] * zc[k];
        }

        float G1 = 1.0f/3.0f;
        float rhof = rhoDebug[index];
        float oneMinusRhof = 1-rhof;
        if(rhof != 0.0f)
        {
            float f1 = 
            R*TT*(
                1.0f + (4.0f*rhof - 2.0f*rhof*rhof)/(oneMinusRhof*oneMinusRhof*oneMinusRhof)
            )-a*rhof - 1.0f/3.0f;
            psxDebug[index] = Mathf.Sqrt(Mathf.Abs(6.0f * rhof * f1 / G1));
            // pressure[index] = rhof/3.0f + G1/6.0f * psx[index]*psx[index]; 
        
            uvwDebug[indexU] /= rhoDebug[index];
            uvwDebug[indexV] /= rhoDebug[index];
            uvwDebug[indexW] /= rhoDebug[index];
        }
        
    }

    void CollisionDebug(int x, int y, int z)
    {
        int k;
        float Tc = 0.3773f*a/(b*R);
        float TT= TT0 *Tc;
        float oneMinusRhow = 1-rhow;
        float f1 = 
        R*TT*(
            1.0f + (4.0f*rhow - 2.0f*rhow*rhow)/(oneMinusRhow*oneMinusRhow*oneMinusRhow)
        )-a*rhow - 1.0f/3.0f;
        float G1 = 1.0f/3.0f;
        float psx_w = Mathf.Sqrt(Mathf.Abs(6.0f * rhow * f1 / G1));

        int index = BufferIndex(x,y,z);
        float[] temp = new float[19];
        int indexU = BufferIndexUvw(index,0);
        int indexV = BufferIndexUvw(index,1);
        int indexW = BufferIndexUvw(index,2);
        float u = uvwDebug[indexU];
        float v = uvwDebug[indexV];
        float w = uvwDebug[indexW];

        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        if(obstDebug[index] != 0)
        {

            for(k = 0; k < 19; k++)
            {
                int indexK = BufferIndexK(index,k);
                temp[k] = ffDebug[indexK];
            }
            for(k = 0; k < 19; k++)
            {
                int indexK = BufferIndexK(index,opp[k]);
                ffDebug[indexK] = temp[k];
            }
        }
        else
        {
            float sumX = 0.0f;
            float sumY = 0.0f;
            float sumZ = 0.0f;
            for(k = 0; k < 19; k++)
            {
                int xTmp = UMod(x + ex[k],lx);
                int yTmp = UMod(y + ey[k],ly);
                int zTmp = UMod(z + ez[k],lz);
                int indexTmp = BufferIndex(xTmp,yTmp,zTmp);
                if(obstDebug[indexTmp] != 0)
                {
                    sumX += t_k[k]*xc[k];
                    sumY += t_k[k]*yc[k];
                    sumZ += t_k[k]*zc[k];
                }
                else
                {
                    Fx += t_k[k]*xc[k] * psxDebug[indexTmp];
                    Fy += t_k[k]*yc[k] * psxDebug[indexTmp];
                    Fz += t_k[k]*zc[k] * psxDebug[indexTmp];
                }
            }

            float Sx = -G1 * sumX * psxDebug[index] * psx_w;
            float Sy = -G1 * sumY * psxDebug[index] * psx_w;
            float Sz = -G1 * sumZ * psxDebug[index] * psx_w;

            Fx = -G1 * psxDebug[index] * Fx;
            Fy = -G1 * psxDebug[index] * Fy;
            Fz = -G1 * psxDebug[index] * Fz;

            u += tauc * (Fx + Sx)/rhoDebug[index];
            v += tauc * (Fy + Sy)/rhoDebug[index];
            w += tauc * (Fz + Sz)/rhoDebug[index];

            float u_squ = u*u + v*v + w*w; 

            for(k = 0; k < 19; k++)
            {
                float u_n = xc[k]*u + yc[k]*v + zc[k]*w;
                float feq = 
                t_k[k]*rhoDebug[index] 
                *(
                    cc*u_n/c_squ + (u_n*cc)*(u_n*cc)/(2.0f*c_squ*c_squ)
                    - u_squ/(2.0f*c_squ)
                ) + t_k[k]*rhoDebug[index];
                int indexK = BufferIndexK(index,k);
                ffDebug[indexK] = feq + (1.0f-1.0f/tauc)*( ffDebug[indexK] - feq );
            }
        }
    }
}
