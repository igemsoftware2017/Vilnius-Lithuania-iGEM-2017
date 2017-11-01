using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;
using Vuforia;

public class TrackingManager : MonoBehaviour, ITrackableEventHandler
{


    public List<GameObject> trackableSubscribeObjects;
    private TrackableBehaviour mTrackableBehaviour;
    private GameObject mainCamera;

    void Start()
    {
        mainCamera = GameObject.FindGameObjectWithTag("GameController");
        mTrackableBehaviour = GetComponent<TrackableBehaviour>();
        if (mTrackableBehaviour)
        {
            mTrackableBehaviour.RegisterTrackableEventHandler(this);
        }
    }

    public void OnTrackableStateChanged(
                                    TrackableBehaviour.Status previousStatus,
                                    TrackableBehaviour.Status newStatus)
    {
        if (newStatus == TrackableBehaviour.Status.DETECTED ||
            newStatus == TrackableBehaviour.Status.TRACKED ||
            newStatus == TrackableBehaviour.Status.EXTENDED_TRACKED)
        {
            //Trigger Target Manager for Detecting Tracking Source Change
            ExecuteEvents.Execute<ITrackableStateHandler>(mainCamera, null,
                    (x, y) => x.OnTrackableFound(gameObject));

            //Trigger Audio Manager and Animations for Stopping Activities
            foreach (var trackableSubscribeObject in trackableSubscribeObjects)
            {
                ExecuteEvents.Execute<ITrackableStateHandler>(trackableSubscribeObject, null, 
                    (x, y) => x.OnTrackableFound(gameObject));
            }           
        }
        else
        {
            //Trigger Audio Manager and Animations for Stopping Activities
            foreach (var trackableSubscribeObject in trackableSubscribeObjects)
            {
                ExecuteEvents.Execute<ITrackableStateHandler>(trackableSubscribeObject, null,
                     (x, y) => x.OnTrackableLost(gameObject));
            }           
        }
    }
 }
