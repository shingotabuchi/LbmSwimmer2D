using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Runtime.InteropServices;

public class LbmSwimmer2D : MonoBehaviour
{
    public Image plotImage;
    public int DIM_X;
    public int DIM_Y;
    public float pr = 0.71f;
    public float ra =   10000.0f;
    public float tauf = 0.8f;
    public float u0 = 0.01f;
    public int loopCount = 1;

    public float minTemp = 0f;
    public float maxTemp = 1f;

    public float minSpeed = 0f;
    public float maxSpeed = 1f;
    public bool debugMode = false;
    public bool debugFrame = false;
    public int particleCount;

    public float particleDensity = 1.25f;
    public float particleRadius;
    public float epsw = 100.0f, zeta = 1.0f;

    Texture2D plotTexture;
    RenderTexture renderTexture;
    int init,collisions,streaming,boundaries,plotSpeed,immersedBoundary,initRoundParticles;
    ComputeBuffer uv,f,force;
    ComputeBuffer roundParticleSmallDataBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterPosBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterVelBuffer;
    ComputeBuffer roundParticleRoundParticlePerimeterFluidVelBuffer;
    ComputeBuffer roundParticleRoundParticleForceOnPerimeterBuffer;
    ComputeBuffer roundParticleInitPosBuffer;
    Vector2[] particleInitPos;
    void OnDestroy()
    {
        uv.Dispose();
        f.Dispose();
        roundParticleSmallDataBuffer.Dispose();
        roundParticleRoundParticlePerimeterPosBuffer.Dispose();
        roundParticleRoundParticlePerimeterVelBuffer.Dispose();
        roundParticleRoundParticlePerimeterFluidVelBuffer.Dispose();
        roundParticleRoundParticleForceOnPerimeterBuffer.Dispose();
    }
    public ComputeShader compute;

    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;
    // Start is called before the first frame update
    void Start()
    {
        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(DIM_X*1080/DIM_Y,1080);
        renderTexture = new RenderTexture(DIM_X,DIM_Y,24);
        renderTexture.enableRandomWrite = true;

        uv = new ComputeBuffer(DIM_X*DIM_Y*2,sizeof(float));
        force = new ComputeBuffer(DIM_X*DIM_Y*2,sizeof(float));
        f = new ComputeBuffer(9*DIM_X*DIM_Y*2,sizeof(float));

        roundParticleSmallDataBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticleSmallData)));
        roundParticleRoundParticlePerimeterPosBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterPos)));
        roundParticleRoundParticlePerimeterVelBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterVel)));
        roundParticleRoundParticlePerimeterFluidVelBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticlePerimeterFluidVel)));
        roundParticleRoundParticleForceOnPerimeterBuffer = new ComputeBuffer(particleCount,Marshal.SizeOf(typeof(RoundParticleForceOnPerimeter)));
        roundParticleInitPosBuffer = new ComputeBuffer(particleCount,sizeof(float) * 2);

        immersedBoundary = compute.FindKernel("ImmersedBoundary");
        compute.SetBuffer(immersedBoundary,"uv",uv);
        compute.SetBuffer(immersedBoundary,"force",force);
        compute.SetBuffer(immersedBoundary,"roundParticleSmallDataBuffer",roundParticleSmallDataBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterPosBuffer",roundParticleRoundParticlePerimeterPosBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterVelBuffer",roundParticleRoundParticlePerimeterVelBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticlePerimeterFluidVelBuffer",roundParticleRoundParticlePerimeterFluidVelBuffer);
        compute.SetBuffer(immersedBoundary,"roundParticleRoundParticleForceOnPerimeterBuffer",roundParticleRoundParticleForceOnPerimeterBuffer);
        compute.SetBuffer(immersedBoundary,"particleInitPos",roundParticleInitPosBuffer);

        compute.SetInt("DIM_X",DIM_X);
        compute.SetInt("DIM_Y",DIM_Y);
        compute.SetFloat("minSpeed",minSpeed);
        compute.SetFloat("maxSpeed",maxSpeed);
        compute.SetFloat("u0",u0);
        compute.SetFloat("tauf",tauf);

        compute.SetInt("particleCount",particleCount);
        compute.SetFloat("particleDensity",particleDensity);
        compute.SetFloat("particleRadius",particleRadius);
        compute.SetFloat("zeta",zeta);
        compute.SetFloat("epsw",epsw);
        particleInitPos = new Vector2[particleCount];
        for(int i = 0; i < particleCount; i++)
        {
            particleInitPos[i] = new Vector2(DIM_X/2,DIM_Y/2);
        }
        roundParticleInitPosBuffer.SetData(particleInitPos);

        initRoundParticles = compute.FindKernel("InitRoundParticles");
        compute.SetBuffer(initRoundParticles,"roundParticleSmallDataBuffer",roundParticleSmallDataBuffer);
        compute.SetBuffer(initRoundParticles,"roundParticleRoundParticlePerimeterPosBuffer",roundParticleRoundParticlePerimeterPosBuffer);
        compute.SetBuffer(initRoundParticles,"roundParticleRoundParticlePerimeterVelBuffer",roundParticleRoundParticlePerimeterVelBuffer);
        compute.SetBuffer(initRoundParticles,"roundParticleRoundParticlePerimeterFluidVelBuffer",roundParticleRoundParticlePerimeterFluidVelBuffer);
        compute.SetBuffer(initRoundParticles,"roundParticleRoundParticleForceOnPerimeterBuffer",roundParticleRoundParticleForceOnPerimeterBuffer);
        compute.SetBuffer(initRoundParticles,"particleInitPos",roundParticleInitPosBuffer);

        init = compute.FindKernel("Init");
        compute.SetBuffer(init,"uv",uv);
        compute.SetBuffer(init,"f",f);

        plotSpeed = compute.FindKernel("PlotSpeed");
        compute.SetBuffer(plotSpeed,"uv",uv);
        compute.SetBuffer(plotSpeed,"f",f);
        compute.SetTexture(plotSpeed,"renderTexture",renderTexture);

        collisions = compute.FindKernel("Collisions");
        compute.SetBuffer(collisions,"f",f);
        compute.SetBuffer(collisions,"uv",uv);
        compute.SetBuffer(collisions,"force",force);

        streaming = compute.FindKernel("Streaming");
        compute.SetBuffer(streaming,"f",f);

        boundaries = compute.FindKernel("Boundaries");
        compute.SetBuffer(boundaries,"f",f);

        compute.Dispatch(init,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        compute.Dispatch(initRoundParticles,(particleCount+63)/64,1,1);
        compute.Dispatch(plotSpeed,(DIM_X+7)/8,(DIM_Y+7)/8,1);

        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
        // compute.Dispatch(plotTemperature,(DIM_X+7)/8,(DIM_Y+7)/8,1);

        // RenderTexture.active = renderTexture;
        // plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        // plotTexture.Apply();
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        int DIM = Mathf.Max(DIM_X,DIM_Y);
        if(debugMode)
        {
            if(debugFrame)
            {
                compute.Dispatch(collisions,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                compute.Dispatch(streaming,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                compute.Dispatch(boundaries,(DIM+63)/64,1,1);
                compute.Dispatch(immersedBoundary,(particleCount+63)/64,1,1);
                debugFrame = false;
            }
        }
        else
        {
            for (int i = 0; i < loopCount; i++)
            {
                compute.Dispatch(collisions,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                compute.Dispatch(streaming,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                compute.Dispatch(boundaries,(DIM+63)/64,1,1);
                compute.Dispatch(immersedBoundary,(particleCount+63)/64,1,1);
            }
        }
        
        compute.Dispatch(plotSpeed,(DIM_X+7)/8,(DIM_Y+7)/8,1);

        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
        
    }

    private void OnValidate() 
    {
        nu = (tauf - 0.5f)/3.0f;
        chi = nu/pr;
        taug = 3.0f*chi + 0.5f;
        rbetag = ra*nu*chi/h/h/h;

        compute.SetFloat("minTemp",minTemp);
        compute.SetFloat("maxTemp",maxTemp);
        compute.SetFloat("minSpeed",minSpeed);
        compute.SetFloat("maxSpeed",maxSpeed);
        compute.SetFloat("u0",u0);
        compute.SetFloat("rbetag",rbetag);
        compute.SetFloat("taug",taug);
        compute.SetFloat("tauf",tauf);
    }
}
