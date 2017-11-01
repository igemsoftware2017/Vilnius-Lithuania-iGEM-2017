using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

public class ScriptMethodExecutor<T> : MonoBehaviour, IAnimationExecutor
{
    private GameObject _baseObject;
    private Action<GameObject, T> _action;
    private T _parameterValue;
    private float _startTime;
    private float? _endTime;

    public ScriptMethodExecutor(GameObject baseObject, Action<GameObject, T> action, T parameterValue,
        float startTime, float? endTime)
    {
        _baseObject = baseObject;
        _action = action;
        _parameterValue = parameterValue;
        _startTime = startTime;
        _endTime = endTime;
    }

    public void Execute()
    {
        _action.Invoke(_baseObject, _parameterValue);
    }

    public float GetStartTime()
    {
        return _startTime;
    }

    public float? GetEndTime()
    {
        return _endTime;
    }

}
