using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Runtime.InteropServices;
using UnityEngine.VFX;
public class Surfactant : MonoBehaviour
{
    public enum PlotMode
    {
        Phase,
        Density,
        Speed,
        Phi,
        C,
    }
    public PlotMode plotMode;
    public Image plotImage;
    public int DIM;
    // public float tauGas=1.0f;
    // public float tauLiq=1.0f;
    public float tau_phi = 1.0f;
    public float tau_rho = 1.0f;
    public int loopCount;
    public float minSpeed = 0f;
    public float maxSpeed = 0.1f;
    public float minRho = 0f;
    public float maxRho = 0.1f;
    public float minPhase = -1f;
    public float maxPhase = 1f;
    public float wallGradient = -0.15f;
    public float sigma = 0.04f;
    public float width = 3f;
    public float initPhaseValue = 1f;
    public float initRho = 1f;
    // public float aconst=0.04f;
    // public float kconst=0.04f;
    public float gammaconst=1.0f;
    public bool debugMode = false;
    public bool debugFrame = false;

    public float aS,bS,gammaS,B0S,B1S,DS;

    private void OnValidate() 
    {
        compute.SetFloat("minSpeed",minSpeed);
        compute.SetFloat("maxSpeed",maxSpeed);
        compute.SetFloat("wallGradient",wallGradient);
        compute.SetFloat("aconst",(3f * sigma*initPhaseValue*initPhaseValue*initPhaseValue*initPhaseValue)/width);
        compute.SetFloat("gammaconst",gammaconst);
        compute.SetFloat("kconst",(3f * sigma * width)/(8f*initPhaseValue*initPhaseValue));
        compute.SetFloat("minRho",minRho);
        compute.SetFloat("maxRho",maxRho);
        compute.SetFloat("minPhase",minPhase);
        compute.SetFloat("maxPhase",maxPhase);
        compute.SetFloat("tau_phi",tau_phi);
        compute.SetFloat("tau_rho",tau_rho);
        compute.SetFloat("initPhase",initPhaseValue);
        compute.SetFloat("initRho",initRho);

        compute.SetFloat("aS",aS);
        compute.SetFloat("bS",bS);
        compute.SetFloat("gammaS",gammaS);
        compute.SetFloat("B0S",B0S);
        compute.SetFloat("B1S",B1S);
        compute.SetFloat("DS",DS);
    }
    Texture2D plotTexture;
    RenderTexture renderTexture;
    int init,wallInit,bulkInit,collisions,streaming,boundaries,plotSpeed,plotPhase;
    int initPhase,plotDensity,plotPhi;
    int cphiStep1,cphiStep2;
    ComputeBuffer uv,f,g,phaseRho,c,phi,muC,muPhi;
    public ComputeShader compute;
    
    private void OnDestroy() {
        uv.Dispose();
        f.Dispose();
        g.Dispose();
        phaseRho.Dispose();
        c.Dispose();
        phi.Dispose();
        muC.Dispose();
        muPhi.Dispose();
    }    

    private void Start() {
        plotTexture = new Texture2D(DIM,DIM);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM,DIM),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(1080,1080);
        
        renderTexture = new RenderTexture(DIM,DIM,24);
        renderTexture.enableRandomWrite = true;

        uv = new ComputeBuffer(DIM*DIM*2,sizeof(float));
        c = new ComputeBuffer(DIM*DIM*2,sizeof(float));
        phi = new ComputeBuffer(DIM*DIM*2,sizeof(float));
        muC = new ComputeBuffer(DIM*DIM*2,sizeof(float));
        muPhi = new ComputeBuffer(DIM*DIM*2,sizeof(float));
        f = new ComputeBuffer(9*DIM*DIM*2,sizeof(float));
        g = new ComputeBuffer(9*DIM*DIM*2,sizeof(float));
        phaseRho = new ComputeBuffer(DIM*DIM*2,sizeof(float));

        compute.SetInt("DIM",DIM);
        compute.SetInt("DIMSqrd",DIM*DIM);
        compute.SetInt("DIMSqrd9",DIM*DIM*9);
        OnValidate();

        init = compute.FindKernel("Init");
        compute.SetBuffer(init,"phaseRho",phaseRho);
        compute.SetBuffer(init,"phi",phi);
        compute.SetBuffer(init,"c",c);

        bulkInit = compute.FindKernel("BulkNodeInit");
        compute.SetBuffer(bulkInit,"uv",uv);
        compute.SetBuffer(bulkInit,"f",f);
        // compute.SetBuffer(bulkInit,"g",g);
        compute.SetBuffer(bulkInit,"phi",phi);
        compute.SetBuffer(bulkInit,"c",c);
        // compute.SetBuffer(bulkInit,"phaseRho",phaseRho);

        wallInit = compute.FindKernel("WallInit");
        compute.SetBuffer(wallInit,"phaseRho",phaseRho);

        initPhase = compute.FindKernel("InitPhase");
        compute.SetBuffer(initPhase,"f",f);
        compute.SetBuffer(initPhase,"g",g);
        compute.SetBuffer(initPhase,"phaseRho",phaseRho);

        plotPhase = compute.FindKernel("PlotPhase");
        compute.SetBuffer(plotPhase,"phaseRho",phaseRho);
        compute.SetTexture(plotPhase,"renderTexture",renderTexture);
        
        plotPhi = compute.FindKernel("PlotPhi");
        compute.SetBuffer(plotPhi,"phi",phi);
        compute.SetTexture(plotPhi,"renderTexture",renderTexture);

        plotSpeed = compute.FindKernel("PlotSpeed");
        compute.SetBuffer(plotSpeed,"uv",uv);
        compute.SetTexture(plotSpeed,"renderTexture",renderTexture);

        plotDensity = compute.FindKernel("PlotDensity");
        compute.SetBuffer(plotDensity,"phaseRho",phaseRho);
        compute.SetTexture(plotDensity,"renderTexture",renderTexture);

        boundaries = compute.FindKernel("BouncebackBoundary");
        compute.SetBuffer(boundaries,"phaseRho",phaseRho);
        compute.SetBuffer(boundaries,"uv",uv);
        compute.SetBuffer(boundaries,"f",f);
        compute.SetBuffer(boundaries,"g",g);

        collisions = compute.FindKernel("Collision");
        // compute.SetBuffer(collisions,"phaseRho",phaseRho);
        compute.SetBuffer(collisions,"phi",phi);
        compute.SetBuffer(collisions,"c",c);
        compute.SetBuffer(collisions,"uv",uv);
        compute.SetBuffer(collisions,"f",f);
        // compute.SetBuffer(collisions,"g",g);

        streaming = compute.FindKernel("Streaming");
        compute.SetBuffer(streaming,"f",f);
        compute.SetBuffer(streaming,"g",g);

        cphiStep1 = compute.FindKernel("CphiStep1");
        compute.SetBuffer(cphiStep1,"uv",uv);
        compute.SetBuffer(cphiStep1,"c",c);
        compute.SetBuffer(cphiStep1,"phi",phi);

        cphiStep2 = compute.FindKernel("CphiStep2");
        compute.SetBuffer(cphiStep2,"c",c);
        compute.SetBuffer(cphiStep2,"phi",phi);

        Init();
        PlotPhase();
    }

    void Init()
    {
        compute.Dispatch(init,(DIM+7)/8,(DIM+7)/8,1);
        // compute.Dispatch(wallInit,(DIM+63)/64,1,1);
        compute.Dispatch(bulkInit,(DIM+7)/8,(DIM+7)/8,1);
    }

    void PlotPhase()
    {
        compute.Dispatch(plotPhase,(DIM+7)/8,(DIM+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
    }
    void PlotDensity()
    {
        compute.Dispatch(plotDensity,(DIM+7)/8,(DIM+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
    }
    void PlotPhi()
    {
        compute.Dispatch(plotPhi,(DIM+7)/8,(DIM+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
    }
    void PlotSpeed()
    {
        compute.Dispatch(plotSpeed,(DIM+7)/8,(DIM+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
    }

    void Step()
    {
        compute.Dispatch(initPhase,(DIM+7)/8,(DIM+7)/8,1);
        // compute.Dispatch(wallInit,(DIM+63)/64,1,1);
        compute.Dispatch(collisions,(DIM+7)/8,(DIM+7)/8,1);
        // compute.Dispatch(boundaries,(DIM+63)/64,1,1);
        // compute.Dispatch(streaming,(DIM+7)/8,(DIM+7)/8,1);
        compute.Dispatch(cphiStep1,(DIM+7)/8,(DIM+7)/8,1);
        compute.Dispatch(cphiStep2,(DIM+7)/8,(DIM+7)/8,1);

    }

    private void FixedUpdate() {
        if(debugMode)
        {
            if(!debugFrame)return;

            debugFrame = false;
        }
        for (int kk = 0; kk < loopCount; kk++)
        {
            Step();
        }
        switch (plotMode)
        {
            case PlotMode.Phase:
            PlotPhase();
            break;
            case PlotMode.Speed:
            PlotSpeed();
            break;
            case PlotMode.Density:
            PlotDensity();
            break;
            case PlotMode.Phi:
            PlotPhi();
            break;
        }
    }
}
