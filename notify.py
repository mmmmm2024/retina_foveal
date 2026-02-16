from plyer import notification

notification.notify(
    title="Pythonで通知",
    message="ここにメッセージを書きます",
    app_name="アプリの名前",
    app_icon="python.ico",
    timeout=10
)