/**
 * Run-completion notifications: a browser notification when a training job
 * finishes while the tab is elsewhere, with a title-prefix fallback where
 * notifications are unavailable or refused. Permission is requested from a
 * user gesture (the Start click), never on load.
 */

export type NotificationState = 'granted' | 'denied' | 'default' | 'unsupported';

const BASE_TITLE = typeof document === 'undefined' ? '' : document.title;

export function notificationState(): NotificationState {
    if (typeof Notification === 'undefined') return 'unsupported';
    return Notification.permission as NotificationState;
}

/** Ask once, from a user gesture; resolves to the resulting state. */
export async function requestNotifications(): Promise<NotificationState> {
    if (typeof Notification === 'undefined') return 'unsupported';
    if (Notification.permission !== 'default') return Notification.permission as NotificationState;
    try {
        return (await Notification.requestPermission()) as NotificationState;
    } catch {
        return notificationState();
    }
}

/** Announce a terminal job. Falls back to a title prefix, which is also what
 * a granted notification gets alongside it while the tab is hidden. */
export function notifyRunFinished(jobId: string, status: string, detail = ''): void {
    const title = `meshio++: run ${jobId} ${status}`;
    if (typeof document !== 'undefined' && document.hidden) {
        document.title = `(${status}) ${BASE_TITLE}`;
        setTimeout(() => {
            if (typeof document !== 'undefined' && !document.hidden) document.title = BASE_TITLE;
        }, 30_000);
    }
    if (typeof Notification === 'undefined' || Notification.permission !== 'granted') return;
    try {
        // eslint-disable-next-line no-new
        new Notification(title, { body: detail, tag: `meshioplusplus-${jobId}` });
    } catch {
        // Some browsers only allow notifications from a service worker; the
        // title fallback above already ran.
    }
}

export function resetTitle(): void {
    if (typeof document !== 'undefined') document.title = BASE_TITLE;
}
