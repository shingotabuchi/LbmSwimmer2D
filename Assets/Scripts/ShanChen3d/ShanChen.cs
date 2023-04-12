using System;
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
    
    int init,plotDensity;
    int stream,calcPressureAndUvw,collision;

    private void OnDestroy() {
        obst.Dispose();
        uvw.Dispose();
        psx.Dispose();
        rho.Dispose();
        ff.Dispose();
    }

    private void Start() 
    {
        plotTexture = new Texture2D(ly,lz);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,ly,lz),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(ly*1080/lz,1080);

        renderTexture = new RenderTexture(ly,lz,24);
        renderTexture.enableRandomWrite = true;

        obst = new ComputeBuffer(lx*ly*lz,sizeof(int));
        uvw = new ComputeBuffer(lx*ly*lz*3,sizeof(float));
        rho = new ComputeBuffer(lx*ly*lz,sizeof(float));
        psx = new ComputeBuffer(lx*ly*lz,sizeof(float));
        ff = new ComputeBuffer(lx*ly*lz*19*2,sizeof(float));

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
        compute.SetFloat("rho_h",rho_h);
        compute.SetFloat("rho_l",rho_l);
        compute.SetFloat("rhow",rhow);
        compute.SetFloat("tauc",tauc);
        compute.SetInt("plotXIndex",plotXIndex);
    }

}
