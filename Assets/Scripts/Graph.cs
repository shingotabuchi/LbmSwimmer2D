using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Graph : MonoBehaviour
{
    public GameObject pointPrefab;
    public Transform plotImage;
    public LineRenderer lineRenderer;
    public Vector2 xRange;
    public Vector2 yRange;
    public bool scale;
    public float slope;
    public Transform stick;
    [Range(0,1)]
    public float useDataPercentage = 1;
    List<Transform> pointList = new List<Transform>();
    public void Plot(float[] xData, float[] yData)
    {
        while(xData.Length * useDataPercentage > pointList.Count)
        {
            pointList.Add(Instantiate(pointPrefab,plotImage).transform);
        }
        while(xData.Length * useDataPercentage < pointList.Count)
        {
            GameObject obj = pointList[pointList.Count-1].gameObject;
            pointList.Remove(pointList[pointList.Count-1]);
            Destroy(obj);
        }

        Vector2 xRange = new Vector2(Mathf.Infinity,-Mathf.Infinity);
        for (int i = 0; i < xData.Length * useDataPercentage; i++)
        {
            xRange.x = Mathf.Min(xData[i],xRange.x);
            xRange.y = Mathf.Max(xData[i],xRange.y);
        }
        Vector2 yRange = new Vector2(Mathf.Infinity,-Mathf.Infinity);
        for (int i = 0; i < xData.Length * useDataPercentage; i++)
        {
            yRange.x = Mathf.Min(yData[i],yRange.x);
            yRange.y = Mathf.Max(yData[i],yRange.y);
        }
        Vector2 originPos = plotImage.GetComponent<RectTransform>().sizeDelta*(-0.5f);
        Vector2 plotSize = plotImage.GetComponent<RectTransform>().sizeDelta;
        for (int i = 0; i < pointList.Count; i++)
        {
            float scaledX = (xData[i] - xRange.x)/(xRange.y - xRange.x);
            float scaledY = (yData[i] - yRange.x)/(yRange.y - yRange.x);
            Vector2 pointPosition = originPos + new Vector2(scaledX * plotSize.x, scaledY * plotSize.y);
            pointPosition = new Vector2(Mathf.Log10(pointPosition.x),Mathf.Log10(pointPosition.y));
            pointList[i].GetComponent<RectTransform>().anchoredPosition = pointPosition;
        }
    }
    public void Plot(double[] xData, double[] yData)
    {
        while(xData.Length * useDataPercentage > pointList.Count)
        {
            pointList.Add(Instantiate(pointPrefab,plotImage).transform);
        }
        while(xData.Length * useDataPercentage < pointList.Count)
        {
            GameObject obj = pointList[pointList.Count-1].gameObject;
            pointList.Remove(pointList[pointList.Count-1]);
            Destroy(obj);
        }

        if(scale)
        {
            xRange = new Vector2(Mathf.Infinity,-Mathf.Infinity);
            for (int i = 1; i < xData.Length * useDataPercentage; i++)
            {
                xRange.x = Mathf.Min(Mathf.Log10((float)xData[i]),xRange.x);
                xRange.y = Mathf.Max(Mathf.Log10((float)xData[i]),xRange.y);
            }
            yRange = new Vector2(Mathf.Infinity,-Mathf.Infinity);
            for (int i = 1; i < xData.Length * useDataPercentage; i++)
            {
                yRange.x = Mathf.Min(Mathf.Log10((float)yData[i]),yRange.x);
                yRange.y = Mathf.Max(Mathf.Log10((float)yData[i]),yRange.y);
            }
        }
        Vector2 originPos = plotImage.GetComponent<RectTransform>().sizeDelta*(-0.5f);
        Vector2 plotSize = plotImage.GetComponent<RectTransform>().sizeDelta;

        lineRenderer.positionCount = pointList.Count-1;
        Vector3[] linePositions = new Vector3[pointList.Count-1];
        for (int i = 1; i < pointList.Count; i++)
        {
            float scaledX = (Mathf.Log10((float)xData[i]) - xRange.x)/(xRange.y - xRange.x);
            float scaledY = (Mathf.Log10((float)yData[i]) - yRange.x)/(yRange.y - yRange.x);
            Vector2 pointPosition = originPos + new Vector2(scaledX * plotSize.x, scaledY * plotSize.y);
            pointList[i].GetComponent<RectTransform>().anchoredPosition = pointPosition;
            linePositions[i-1] = pointList[i].position;
        }
        lineRenderer.SetPositions(linePositions);

        stick.transform.rotation = Quaternion.Euler(0,0,Mathf.Atan(slope /((yRange.y - yRange.x)/(xRange.y - xRange.x))) * Mathf.Rad2Deg);
    }
}
