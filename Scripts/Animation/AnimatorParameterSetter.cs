using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AnimatorParameterSetter<T> : MonoBehaviour, IAnimationExecutor {

    private Animator _animator;
    private string _parameterName;
    private T _parameterValue;
    private float _startTime;
    private float? _endTime;

    public AnimatorParameterSetter(Animator animator, string parameterName, T parameterValue,
        float startTime, float? endTime)
    {
        _animator = animator;
        _parameterName = parameterName;
        _parameterValue = parameterValue;
        _startTime = startTime;
        _endTime = endTime;
    }

    public void Execute()
    {
        if(typeof(float) == typeof(T))
        {
            _animator.SetFloat(_parameterName, (float)(object)_parameterValue);
        }else if (typeof(int) == typeof(T))
        {
            _animator.SetInteger(_parameterName, (int)(object)_parameterValue);
        }
        else if (typeof(bool) == typeof(T))
        {
            _animator.SetBool(_parameterName, (bool)(object)_parameterValue);
        }
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
