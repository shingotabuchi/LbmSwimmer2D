using System.IO;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SwimmerDataWriter : MonoBehaviour
{
    [SerializeField] LbmSwimmer2D simulation;
    [SerializeField] int[] timeSteps;
    [SerializeField] bool distandtime;
    [SerializeField] int totalTime;
    int dataIndex = 0;
    string timeanddist = "";
    string precision = "F10";
    void Update()
    {
        if(!distandtime)
        {
            if(dataIndex>=timeSteps.Length)
            {
                UnityEditor.EditorApplication.isPlaying = false;
                return;
            }
            if(simulation.timeFrame >= timeSteps[dataIndex])
            {
                simulation.FillUvBuffer();
                simulation.FillParticleBuffer();
                WriteField();
                dataIndex++;
            }
        }
        else
        {
            if(dataIndex>=totalTime)
            {
                WriteDistTime();
                UnityEditor.EditorApplication.isPlaying = false;
                return;
            }
            if(simulation.timeFrame >= dataIndex)
            {
                simulation.FillParticleBuffer();
                // Vector2[] offsets = new Vector2[4]{new Vector2(simulation.DIM_X,0),new Vector2(-simulation.DIM_X,0),new Vector2(0,simulation.DIM_Y),new Vector2(0,-simulation.DIM_Y)};
                float dist = (simulation.debugSmallData[0].pos - simulation.debugSmallData[1].pos).magnitude;
                // for (int i = 0; i < 4; i++)
                // {
                //     dist = Mathf.Min(dist,(simulation.debugSmallData[0].pos - simulation.debugSmallData[1].pos + offsets[i]).magnitude);
                // }                
                // timeanddist += simulation.timeFrame.ToString() + " "+ dist.ToString(precision) + "\n";
                timeanddist += dist.ToString(precision) + "\n";
                dataIndex++;
            }
        }
        
    }
    void WriteDistTime()
    {
        float beta = simulation.squirmerBeta;
        int time = simulation.timeFrame;
        string parameters = beta.ToString("0") + "_" + time.ToString();
        DateTime today = DateTime.Today;
        string path = @"Assets/Scripts/python/" + today.ToString("yyyy/mm/dd") + "/";
        Directory.CreateDirectory(path);

        File.WriteAllText(path + "timeanddist" + parameters + ".txt", timeanddist);
    }
    void WriteField()
    {
        string U = "";
        string V = "";
        string px = "";
        string py= "";
        string theta= "";
        string velx= "";
        string vely= "";
        
        for (int j = 0; j < simulation.DIM_Y; j++)
        {
            for (int i = 0; i < simulation.DIM_X; i++)
            {
                int index = i + j*simulation.DIM_X;
                U += simulation.uvBuffer[index*2 + 0].ToString(precision);
                V += simulation.uvBuffer[index*2 + 1].ToString(precision);

                if(i < simulation.DIM_X-1)
                {
                    U += " ";
                    V += " ";
                }
            }
            if(j < simulation.DIM_Y-1)
            {
                U += "\n";
                V += "\n";
            }
        }
        for (int i = 0; i < simulation.particleCount; i++)
        {
            px += simulation.debugSmallData[i].pos.x.ToString(precision);
            py += simulation.debugSmallData[i].pos.y.ToString(precision);
            velx += simulation.debugSmallData[i].vel.x.ToString(precision);
            vely += simulation.debugSmallData[i].vel.y.ToString(precision);
            theta += simulation.debugSmallData[i].theta.ToString(precision);
            if(i < simulation.particleCount-1)
            {
                px += "\n";
                py += "\n";
                velx += "\n";
                vely += "\n";
                theta += "\n";
            }
        }
        float beta = simulation.squirmerBeta;
        int time = simulation.timeFrame;
        if(time == 1) time = 0;
        string parameters = beta.ToString("0") + "_" + time.ToString();
        DateTime today = DateTime.Today;
        string path = @"Assets/Scripts/python/" + today.ToString("yyyy/mm/dd") + "/";
        Directory.CreateDirectory(path);

        File.WriteAllText(path + "U" + parameters + ".txt", U);
        File.WriteAllText(path + "V" + parameters + ".txt", V);
        File.WriteAllText(path + "px" + parameters + ".txt", px);
        File.WriteAllText(path + "py" + parameters + ".txt", py);
        File.WriteAllText(path + "velx" + parameters + ".txt", velx);
        File.WriteAllText(path + "vely" + parameters + ".txt", vely);
        File.WriteAllText(path + "theta" + parameters + ".txt", theta);
    }
}
