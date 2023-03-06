using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.Runtime.InteropServices;
using UnityEngine.VFX;
using FftSharp;
public class LbmNoSwimmer2D : MonoBehaviour
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

    Texture2D plotTexture;
    RenderTexture renderTexture;
    int init,collisions,streaming,boundaries,plotSpeed,immersedBoundary;
    ComputeBuffer uv,f,force;
    Vector2[] debugFluidVel;

    public VisualEffect vfx;
    GraphicsBuffer velocityGraphicsBuffer;
    public bool showGraph;
    public Graph graph;
    float[] uvBuffer;

    bool isTouched;
    [Range(0,1)]
    public float touchForceRadius;
    public float touchForce;


    void OnDestroy()
    {
        uv.Dispose();
        f.Dispose();
        force.Dispose();
        velocityGraphicsBuffer.Dispose();
    }
    public ComputeShader compute;

    float umax, umin, tmp, u2, nu, chi, norm, taug, rbetag, h;

    double[] qArray;
    Dictionary<int,int> gridPosToQIndex = new Dictionary<int,int>();
    // Start is called before the first frame update
    void Start()
    {
        if(showGraph) graph.transform.gameObject.SetActive(true);
        else graph.transform.gameObject.SetActive(false);
        
        double[] qFromX = FftSharp.Transform.FFTfreq(1, DIM_X,false);
        double[] qFromY = FftSharp.Transform.FFTfreq(1, DIM_Y,false);
        HashSet<double> qSet = new HashSet<double>();
        List<double> qList = new List<double>();
        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                double q = Math.Sqrt(qFromX[i]*qFromX[i] + qFromY[j]*qFromY[j]);
                if(!qSet.Contains(q))
                {
                    qList.Add(q);
                }   
                qSet.Add(q);
            }
        }

        qList.Sort();
        for (int i = 0; i < DIM_X; i++)
        {
            for (int j = 0; j < DIM_Y; j++)
            {
                double q = Math.Sqrt(qFromX[i]*qFromX[i] + qFromY[j]*qFromY[j]);
                for (int k = 0; k < qList.Count; k++)
                {
                    if(q == qList[k])
                    {
                        gridPosToQIndex.Add(i + j*DIM_X,k);
                        break;
                    }
                }
            }
        }

        qArray = qList.ToArray();

        velocityGraphicsBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Structured, DIM_X*DIM_Y, sizeof(float)*2);
        vfx.SetGraphicsBuffer("VelocityBuffer",velocityGraphicsBuffer);
        vfx.SetInt("DIM_X",DIM_X);
        vfx.SetInt("DIM_Y",DIM_Y);


        plotTexture = new Texture2D(DIM_X,DIM_Y);
        plotTexture.filterMode = FilterMode.Point;
        plotImage.sprite = Sprite.Create(plotTexture, new Rect(0,0,DIM_X,DIM_Y),UnityEngine.Vector2.zero);
        ((RectTransform)plotImage.transform).sizeDelta = new Vector2(DIM_X*1080/DIM_Y,1080);

        renderTexture = new RenderTexture(DIM_X,DIM_Y,24);
        renderTexture.enableRandomWrite = true;

        uv = new ComputeBuffer(DIM_X*DIM_Y*2,sizeof(float));
        force = new ComputeBuffer(DIM_X*DIM_Y*2,sizeof(float));
        f = new ComputeBuffer(9*DIM_X*DIM_Y*2,sizeof(float));
        uvBuffer = new float[DIM_X*DIM_Y*2];

        immersedBoundary = compute.FindKernel("ImmersedBoundary");
        compute.SetBuffer(immersedBoundary,"uv",uv);
        compute.SetBuffer(immersedBoundary,"force",force);

        compute.SetInt("DIM_X",DIM_X);
        compute.SetInt("DIM_Y",DIM_Y);
        compute.SetFloat("minSpeed",minSpeed);
        compute.SetFloat("maxSpeed",maxSpeed);
        compute.SetFloat("u0",u0);
        compute.SetFloat("touchForce",touchForce);
        compute.SetFloat("touchForceRadius",touchForceRadius);
        compute.SetFloat("tauf",tauf);

        init = compute.FindKernel("Init");
        compute.SetBuffer(init,"force",force);
        compute.SetBuffer(init,"uv",uv);
        compute.SetBuffer(init,"f",f);

        plotSpeed = compute.FindKernel("PlotSpeed");
        compute.SetBuffer(plotSpeed,"uv",uv);
        compute.SetBuffer(plotSpeed,"f",f);
        compute.SetTexture(plotSpeed,"renderTexture",renderTexture);

        collisions = compute.FindKernel("Collisions");
        compute.SetBuffer(collisions,"velocityBuffer",velocityGraphicsBuffer);
        compute.SetBuffer(collisions,"f",f);
        compute.SetBuffer(collisions,"uv",uv);
        compute.SetBuffer(collisions,"force",force);

        streaming = compute.FindKernel("Streaming");
        compute.SetBuffer(streaming,"f",f);

        boundaries = compute.FindKernel("Boundaries");
        compute.SetBuffer(boundaries,"f",f);

        compute.Dispatch(init,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        compute.Dispatch(plotSpeed,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
        // compute.Dispatch(plotTemperature,(DIM_X+7)/8,(DIM_Y+7)/8,1);

        // RenderTexture.active = renderTexture;
        // plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        // plotTexture.Apply();
    }

    void Update()
    {
        if(Input.GetMouseButton(0))
        {
            Vector2 touchPos = ((Vector2)Input.mousePosition - new Vector2(1920f,1080f)/2f - GetComponent<RectTransform>().anchoredPosition)/1080f;
            if(Mathf.Abs(touchPos.x) <= 1f && Mathf.Abs(touchPos.y) <= 1f)
            {
                if(!isTouched)
                {
                    touchForce = Mathf.Abs(touchForce);
                    compute.SetFloat("touchForce",touchForce);
                    isTouched = true;
                    compute.SetBool("isTouched",isTouched);
                }
                compute.SetVector("touchTextureCoord",touchPos);
            }
            else if(isTouched)
            {
                isTouched = false;
                compute.SetBool("isTouched",isTouched);
            }
        }
        else if(Input.GetMouseButton(1))
        {
            Vector2 touchPos = ((Vector2)Input.mousePosition - new Vector2(1920f,1080f)/2f - GetComponent<RectTransform>().anchoredPosition)/1080f;
            if(Mathf.Abs(touchPos.x) <= 1f && Mathf.Abs(touchPos.y) <= 1f)
            {
                if(!isTouched)
                {
                    touchForce = -Mathf.Abs(touchForce);
                    compute.SetFloat("touchForce",touchForce);
                    isTouched = true;
                    compute.SetBool("isTouched",isTouched);
                }
                compute.SetVector("touchTextureCoord",touchPos);
            }
            else if(isTouched)
            {
                isTouched = false;
                compute.SetBool("isTouched",isTouched);
            }
        }
        else if(isTouched)
        {
            isTouched = false;
            compute.SetBool("isTouched",isTouched);
        }
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
                // compute.Dispatch(boundaries,(DIM+63)/64,1,1);
                foreach(Vector2 vel in debugFluidVel)
                {
                    print(vel);
                }
                debugFrame = false;
            }
        }
        else
        {
            for (int i = 0; i < loopCount; i++)
            {
                compute.Dispatch(collisions,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                compute.Dispatch(streaming,(DIM_X+7)/8,(DIM_Y+7)/8,1);
                // compute.Dispatch(boundaries,(DIM+63)/64,1,1);
            }
        }
        vfx.SetFloat("VelocityScale",(10f/(float)DIM_Y) * loopCount / Time.deltaTime);
        compute.Dispatch(plotSpeed,(DIM_X+7)/8,(DIM_Y+7)/8,1);
        RenderTexture.active = renderTexture;
        plotTexture.ReadPixels(new Rect(0, 0, renderTexture.width, renderTexture.height), 0, 0);
        plotTexture.Apply();
        
        if(showGraph)
        {
            uv.GetData(uvBuffer);
            // row fft
            List<Complex[]> uColumnsFfted = new List<Complex[]>();
            List<Complex[]> vColumnsFfted = new List<Complex[]>();
            for (int j = 0; j < DIM_Y; j++)
            {
                double[] uRow = new double[DIM_X];
                double[] vRow = new double[DIM_X];
                for (int i = 0; i < DIM_X; i++)
                {
                    uRow[i] = (double)uvBuffer[(i + j*DIM_X) * 2 + 0];
                    vRow[i] = (double)uvBuffer[(i + j*DIM_X) * 2 + 1];
                }
                Complex[] uRowFfted = FftSharp.Transform.FFT(uRow);
                Complex[] vRowFfted = FftSharp.Transform.FFT(vRow);
                for (int i = 0; i < DIM_X; i++)
                {
                    if(j==0)
                    {
                        uColumnsFfted.Add(new Complex[DIM_Y]);
                        vColumnsFfted.Add(new Complex[DIM_Y]);
                    }
                    uColumnsFfted[i][j] = uRowFfted[i];
                    vColumnsFfted[i][j] = vRowFfted[i];
                }
                if(uRowFfted.Length != DIM_X) print("wtf");
            }
            double[] energySpectrum = new double[qArray.Length];
            for (int i = 0; i < energySpectrum.Length; i++)
            {
                energySpectrum[i] = 0.0;
            }
            for (int i = 0; i < DIM_X; i++)
            {
                FftSharp.Transform.FFT(uColumnsFfted[i]);
                FftSharp.Transform.FFT(vColumnsFfted[i]);
                for (int j = 0; j < DIM_Y; j++)
                {
                    energySpectrum[gridPosToQIndex[i + j*DIM_X]] += uColumnsFfted[i][j].MagnitudeSquared + vColumnsFfted[i][j].MagnitudeSquared;
                }
            }

            graph.Plot(qArray,energySpectrum);
        }
    }

    private void OnValidate() 
    {
        compute.SetFloat("minSpeed",minSpeed);
        compute.SetFloat("maxSpeed",maxSpeed);
        compute.SetFloat("u0",u0);
        compute.SetFloat("tauf",tauf);
        compute.SetFloat("touchForce",touchForce);
        compute.SetFloat("touchForceRadius",touchForceRadius);

        if(showGraph) graph.transform.gameObject.SetActive(true);
        else graph.transform.gameObject.SetActive(false);
    }
}
