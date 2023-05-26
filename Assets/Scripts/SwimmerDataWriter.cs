using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SwimmerDataWriter : MonoBehaviour
{
    [SerializeField] LbmSwimmer2D simulation;
    [SerializeField] int[] timeSteps;
    int dataIndex = 0;
    void Update()
    {
        if(dataIndex>=timeSteps.Length) return;
        if(simulation.timeFrame >= timeSteps[dataIndex])
        {
            simulation.FillUvBuffer();
            WriteField();
            dataIndex++;
        }
    }

    void WriteField()
    {
        string U = "";
        string V = "";
        string px = "";
        string py= "";
        for (int j = 0; j < simulation.DIM_Y; j++)
        {
            for (int i = 0; i < simulation.DIM_X; i++)
            {
                int index = i + j*simulation.DIM_X;
                U += simulation.uvBuffer[index*2 + 0].ToString("0.0000000000000000000");
                V += simulation.uvBuffer[index*2 + 1].ToString("0.0000000000000000000");

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
            px += simulation.debugSmallData[i].pos.x.ToString("0.0000000000000000000");
            py += simulation.debugSmallData[i].pos.y.ToString("0.0000000000000000000");
            if(i < simulation.particleCount-1)
            {
                px += "\n";
                py += "\n";
            }
        }
        float beta = simulation.squirmerBeta*2f;
        int time = simulation.timeFrame;
        if(time == 1) time = 0;
        string parameters = beta.ToString("0") + "_" + time.ToString();
        File.WriteAllText(@"Assets/Scripts/python/U" + parameters + ".txt", U);
        File.WriteAllText(@"Assets/Scripts/python/V" + parameters + ".txt", V);
        File.WriteAllText(@"Assets/Scripts/python/px" + parameters + ".txt", px);
        File.WriteAllText(@"Assets/Scripts/python/py" + parameters + ".txt", py);
    }
}
